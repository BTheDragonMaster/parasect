import React, { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import Graph from 'graphology';
import Sigma from 'sigma';
import forceAtlas2 from 'graphology-layout-forceatlas2';
import {
    Box, Typography, Slider, FormControl, InputLabel, Select, MenuItem, Autocomplete,
    TextField, List, ListItemButton, ListItemText, Chip, Button, CircularProgress, Paper, Divider,
    FormControlLabel, Switch, Alert, IconButton, Tooltip, createFilterOptions, useMediaQuery, useTheme,
} from '@mui/material';
import DownloadIcon from '@mui/icons-material/Download';
import KeyboardDoubleArrowLeftIcon from '@mui/icons-material/KeyboardDoubleArrowLeft';
import TuneIcon from '@mui/icons-material/Tune';

import Loading from '../components/Loading';
import { DistanceBadge, DistanceLegend } from '../components/DistanceBadge';
import { useColorMode } from '../theme/ColorModeContext';
import { categoricalColor, MUTED_MARK, OTHER, SAFE_SLOTS } from '../theme';
import { graphToSvg, svgToPng } from '../utils/networkExport';
import { createZip, downloadBlob } from '../utils/zip';

/** Colour modes backed by a database-wide category list (/api/network/categories).
 * 'cluster' is not: its categories only exist for the current threshold, so they
 * are derived from the graph instead. */
const CATEGORY_FIELDS = ['substrate', 'genus', 'kingdom'];

/** Matches the number of validated, colour-blind-safe categorical slots. */
const DEFAULT_LEGEND_SIZE = SAFE_SLOTS;

/** MUI's Autocomplete mounts every matching option, so the listbox is capped:
 * colouring by cluster at a low threshold offers thousands of them, and a
 * growing reference database will do the same to the substrate list. The lists
 * are frequency-ordered, so the head is the useful part, and the tail is reported
 * rather than silently dropped. */
const OPTION_RENDER_LIMIT = 100;
const MORE_OPTIONS = '\u0000more';
const filterCategories = createFilterOptions({ trim: true });

/** Sidebar width bounds. The floor is what the legend chips and the threshold
 * slider need to stay usable; the ceiling stops a drag from squeezing the graph
 * out of its own page. */
const SIDEBAR_MIN = 280;
const SIDEBAR_MAX = 640;
const SIDEBAR_DEFAULT = 340;
/** How far an arrow key moves the divider, for resizing without a pointer. */
const SIDEBAR_STEP = 24;
const SIDEBAR_STORAGE_KEY = 'paras-network-sidebar';

const clampSidebar = (width) => Math.min(SIDEBAR_MAX, Math.max(SIDEBAR_MIN, Math.round(width)));

/** Domain names are one long token with no spaces, like "Q04747.3.A2|Q45675.A2|
 * WP_069322367.1.A2", so the default word-break leaves them hanging out of
 * whatever column they are in. They have to be allowed to break mid-token. */
const WRAP_ANYWHERE = { overflowWrap: 'anywhere' };

/** Sidebar width and open state from the last visit, or the defaults. */
function readSidebarPreference() {
    try {
        const stored = JSON.parse(window.localStorage.getItem(SIDEBAR_STORAGE_KEY) || '{}');
        return {
            width: Number.isFinite(stored.width) ? clampSidebar(stored.width) : SIDEBAR_DEFAULT,
            open: stored.open !== false,
        };
    } catch (err) {
        // blocked or corrupted storage: the defaults are fine, and the next
        // resize will try to write again
        return { width: SIDEBAR_DEFAULT, open: true };
    }
}

/**
 * Give each legend entry a colour slot, reusing freed slots for new entries.
 *
 * Colour has to follow the category, not its position in the list: if removing
 * "leucine" shifted every colour after it, the whole picture would repaint and
 * nothing you had learned from it would still hold. Entries therefore keep the
 * slot they were given until they leave the legend, and a new entry takes the
 * lowest slot nobody is using.
 *
 * @param {Map<string, number>} current - existing label -> slot assignments.
 * @param {string[]} labels - the labels that should have a slot.
 * @returns {Map<string, number>} - assignments covering exactly `labels`.
 */
function assignSlots(current, labels) {
    const next = new Map();
    const taken = new Set();
    labels.forEach((label) => {
        if (current.has(label)) {
            next.set(label, current.get(label));
            taken.add(current.get(label));
        }
    });
    labels.forEach((label) => {
        if (next.has(label)) return;
        let slot = 0;
        while (taken.has(slot)) slot += 1;
        taken.add(slot);
        next.set(label, slot);
    });
    return next;
}

/**
 * Hover marker for a node: a ring, and no text.
 *
 * sigma's stock hover renderer paints a hardcoded white pill and repeats the
 * node's label inside it which turns into a second, unreadable copy of the name sitting on
 * top of the tooltip that already shows it (and a white box on a dark page).
 * Drawing just a ring keeps leaves the naming to the tooltip.
 *
 * @param {CanvasRenderingContext2D} context - the hover layer.
 * @param {object} data - reduced node display data.
 * @param {object} settings - sigma settings; labelColor tracks the colour mode.
 * @returns {void}
 */
function drawNodeHoverRing(context, data, settings) {
    context.strokeStyle = settings.labelColor?.color || '#000';
    context.lineWidth = 2;
    context.beginPath();
    context.arc(data.x, data.y, data.size + 3, 0, Math.PI * 2);
    context.closePath();
    context.stroke();
}

/** A single domain's labels as the same {label: domainCount} shape the server sends
 * for clusters, so one code path can count and match both kinds of node. */
function compositionOf(labels) {
    const list = (labels || []).filter(Boolean);
    if (!list.length) return { unknown: 1 };
    return list.reduce((acc, label) => ({ ...acc, [label]: 1 }), {});
}

/** The highlighted label this node carries most domains of, or null if it carries none. */
function bestMatch(composition, highlighted) {
    let best = null;
    let bestCount = 0;
    for (const label of highlighted) {
        const count = composition[label] || 0;
        if (count > bestCount) {
            best = label;
            bestCount = count;
        }
    }
    return best;
}

/** "leucine 12, valine 3, +4 more"; a cluster's composition, biggest share first. */
function describeComposition(composition, limit = 4) {
    const entries = Object.entries(composition || {}).sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0]));
    if (!entries.length) return 'unknown';
    const shown = entries.slice(0, limit).map(([label, count]) => `${label} ${count}`).join(', ');
    return entries.length > limit ? `${shown}, +${entries.length - limit} more` : shown;
}

/**
 * Sequence-similarity network of the reference PARAS/PARASECT database.
 *
 * Domains are grouped into clusters (connected components at a Hamming-distance
 * threshold over their 34-residue extended signatures) to keep the graph
 * renderable. Click a cluster to expand it into its individual domains. A
 * cluster is coloured by the majority vote of its members, and every node also
 * carries its full composition so the legend can highlight a category that is
 * present without winning that vote.
 *
 * @returns {React.ReactElement} - The network graph page.
 */
const NetworkGraph = () => {
    const containerRef = useRef(null);
    const sigmaRef = useRef(null);
    const graphRef = useRef(null);
    const sidebarRef = useRef(null);
    const theme = useTheme();
    const { mode: colorMode } = useColorMode();
    // below md the sidebar stacks above the graph at full width, where a
    // vertical divider has nothing to drag
    const canResizeSidebar = useMediaQuery(theme.breakpoints.up('md'));

    // refs mirroring state that the sigma reducers (set up once) need to read live
    const colorByRef = useRef('substrate');
    const highlightRef = useRef([]);
    const searchHighlightRef = useRef(null);
    const paintRef = useRef(() => '#999999');
    // the search marker is a reserved status colour, never a categorical slot
    const searchMarkRef = useRef('#B3261E');
    const edgeColorRef = useRef('#dddddd');
    const themeRef = useRef(null);

    const [displayThreshold, setDisplayThreshold] = useState(5);
    const [threshold, setThreshold] = useState(5);
    const [colorBy, setColorBy] = useState('substrate');
    const [loading, setLoading] = useState(true);
    const [meta, setMeta] = useState({ total_domains: 0, cluster_count: 0 });
    const [hoverInfo, setHoverInfo] = useState(null);

    // every category in the database per colour field, most frequent first
    const [categories, setCategories] = useState({});
    // 'cluster' mode equivalent, derived from whatever is currently in the graph
    const [clusterCategories, setClusterCategories] = useState([]);
    // which categories are in the legend, per colour mode; undefined = not seeded yet
    const [legend, setLegend] = useState({});
    // legend entries the user clicked: only nodes containing one of these stay coloured
    const [highlighted, setHighlighted] = useState([]);
    const [matchStats, setMatchStats] = useState(null);
    // label -> colour slot, per colour mode; see assignSlots
    const [slots, setSlots] = useState({});
    const [showLabels, setShowLabels] = useState(true);
    const [exporting, setExporting] = useState(false);

    const [searchQuery, setSearchQuery] = useState('');
    const [searchMatches, setSearchMatches] = useState([]);
    const [neighborResult, setNeighborResult] = useState(null);
    const [searchLoading, setSearchLoading] = useState(false);

    const [sidebar, setSidebar] = useState(readSidebarPreference);

    // the graph is the page; a width the user dragged to, or a sidebar they
    // closed to see it, should still be there next time
    useEffect(() => {
        try {
            window.localStorage.setItem(SIDEBAR_STORAGE_KEY, JSON.stringify(sidebar));
        } catch (err) {
            // the layout still applies for this session; nothing to recover
        }
    }, [sidebar]);

    const setSidebarWidth = useCallback(
        (width) => setSidebar((prev) => ({ ...prev, width: clampSidebar(width) })),
        [],
    );

    /**
     * Drag the divider between the sidebar and the graph.
     *
     * The listeners go on the window rather than the divider so the drag
     * survives the pointer outrunning a 6px target, and `setPointerCapture`
     * would not help: the canvas swallows the events it captures. Text
     * selection is suppressed for the duration, or dragging left selects the
     * whole sidebar.
     */
    const startSidebarResize = useCallback((event) => {
        event.preventDefault();
        const startX = event.clientX;
        const startWidth = sidebarRef.current?.getBoundingClientRect().width || SIDEBAR_DEFAULT;
        const onMove = (move) => setSidebarWidth(startWidth + move.clientX - startX);
        const onUp = () => {
            window.removeEventListener('pointermove', onMove);
            window.removeEventListener('pointerup', onUp);
            document.body.style.removeProperty('cursor');
            document.body.style.removeProperty('user-select');
        };
        document.body.style.cursor = 'col-resize';
        document.body.style.userSelect = 'none';
        window.addEventListener('pointermove', onMove);
        window.addEventListener('pointerup', onUp);
    }, [setSidebarWidth]);

    const onSidebarResizeKey = useCallback((event) => {
        if (event.key !== 'ArrowLeft' && event.key !== 'ArrowRight') return;
        event.preventDefault();
        setSidebarWidth(
            (sidebarRef.current?.getBoundingClientRect().width || SIDEBAR_DEFAULT)
            + (event.key === 'ArrowRight' ? SIDEBAR_STEP : -SIDEBAR_STEP),
        );
    }, [setSidebarWidth]);

    // memoised so the empty-list fallback doesn't hand every dependent hook a
    // fresh array identity on each render
    const availableCategories = useMemo(
        () => (colorBy === 'cluster' ? clusterCategories : (categories[colorBy] || [])),
        [colorBy, clusterCategories, categories],
    );
    // memoised for the same reason as availableCategories: the `|| []` fallback
    // would otherwise be a new array on every render
    const legendLabels = useMemo(() => legend[colorBy] || [], [legend, colorBy]);

    const countByLabel = useMemo(
        () => new Map(availableCategories.map((c) => [c.label, c.count])),
        [availableCategories],
    );

    // keep a slot for every legend entry, and only for legend entries
    useEffect(() => {
        setSlots((prev) => {
            const current = prev[colorBy] || new Map();
            const next = assignSlots(current, legendLabels);
            const unchanged = next.size === current.size
                && [...next].every(([label, slot]) => current.get(label) === slot);
            return unchanged ? prev : { ...prev, [colorBy]: next };
        });
    }, [colorBy, legendLabels]);

    const activeSlots = useMemo(() => slots[colorBy] || new Map(), [slots, colorBy]);

    /** A legend entry's colour, or null when the category isn't in the legend. */
    const colorForLabel = useCallback((label) => {
        const slot = activeSlots.get(label);
        return slot === undefined ? null : categoricalColor(slot, colorMode);
    }, [activeSlots, colorMode]);

    const otherColor = OTHER[colorMode];
    const mutedColor = MUTED_MARK[colorMode];

    /**
     * The colour a node takes, shared by the renderer and the SVG export so the
     * downloaded figure is the picture on screen.
     */
    const paintNode = useCallback((node, attrs) => {
        const mode = colorByRef.current;
        const activeHighlights = highlightRef.current;
        if (activeHighlights.length) {
            const match = bestMatch(attrs.composition[mode] || {}, activeHighlights);
            return match === null ? mutedColor : (colorForLabel(match) || otherColor);
        }
        return colorForLabel(attrs.dominant[mode] || 'unknown') || otherColor;
    }, [colorForLabel, mutedColor, otherColor]);

    useEffect(() => { paintRef.current = paintNode; sigmaRef.current?.refresh(); }, [paintNode]);

    const edgeColor = theme.palette.surface.border;
    const surfaceColor = theme.palette.background.paper;

    // written during render so the one-time sigma setup below, which runs after
    // this render but before any effect that depends on `theme`, sees the real one
    themeRef.current = theme;
    edgeColorRef.current = edgeColor;

    // sigma keeps its own copy of these, so push them when the mode flips
    useEffect(() => {
        searchMarkRef.current = theme.palette.error.main;
        edgeColorRef.current = edgeColor;
        const renderer = sigmaRef.current;
        const g = graphRef.current;
        if (!renderer || !g) return;
        renderer.setSetting('labelColor', { color: theme.palette.text.primary });
        renderer.setSetting('edgeColor', 'default');
        renderer.setSetting('defaultEdgeColor', edgeColor);
        g.forEachEdge((edge) => g.setEdgeAttribute(edge, 'color', edgeColor));
        renderer.refresh();
    }, [theme, edgeColor]);

    // labels off makes for a much cleaner figure; sigma renders them itself
    useEffect(() => {
        sigmaRef.current?.setSetting('renderLabels', showLabels);
        sigmaRef.current?.refresh();
    }, [showLabels]);

    /** Cluster categories and highlight match counts both need a pass over the graph. */
    const recomputeFromGraph = useCallback(() => {
        const g = graphRef.current;
        if (!g) return;

        const clusterSizes = new Map();
        let matchedNodes = 0;
        let matchedDomains = 0;
        const activeHighlights = highlightRef.current;
        const mode = colorByRef.current;

        g.forEachNode((node, attrs) => {
            Object.entries(attrs.composition.cluster).forEach(([key, count]) => {
                clusterSizes.set(key, (clusterSizes.get(key) || 0) + count);
            });
            if (activeHighlights.length) {
                const composition = attrs.composition[mode] || {};
                const domains = activeHighlights.reduce((sum, label) => sum + (composition[label] || 0), 0);
                if (domains > 0) {
                    matchedNodes += 1;
                    matchedDomains += domains;
                }
            }
        });

        setClusterCategories(
            [...clusterSizes.entries()]
                .sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0]))
                .map(([label, count]) => ({ label, count })),
        );
        setMatchStats(
            activeHighlights.length
                ? { nodes: matchedNodes, domains: matchedDomains, total: g.order }
                : null,
        );
    }, []);

    const refreshSigma = useCallback(() => {
        recomputeFromGraph();
        sigmaRef.current?.refresh();
    }, [recomputeFromGraph]);

    // one-time sigma setup
    useEffect(() => {
        const graph = new Graph();
        graphRef.current = graph;

        const renderer = new Sigma(graph, containerRef.current, {
            minCameraRatio: 0.02,
            maxCameraRatio: 3,
            labelRenderedSizeThreshold: 6,
            // sigma defaults labels to #000, and the effect that themes them runs
            // before this one on mount, so the mount-time theme has to be passed
            // here or a first load in dark mode draws black labels on a black page
            labelColor: { color: themeRef.current.palette.text.primary },
            defaultEdgeColor: edgeColorRef.current,
            defaultDrawNodeHover: drawNodeHoverRing,
            // the container is still being laid out on the first paint in a narrow
            // window, and sigma throws on a zero-width container rather than waiting
            // for one; an uncaught error there takes the whole page down
            allowInvalidContainer: true,
        });
        sigmaRef.current = renderer;

        // sigma only watches the window for resizes, so it never notices the
        // container growing on its own: not when it starts at zero width above, and
        // not when the sidebar wraps or a split pane is dragged. refresh() resizes
        // the canvases itself, so handing it every container resize covers both.
        const resizeObserver = new ResizeObserver(() => renderer.scheduleRefresh());
        resizeObserver.observe(containerRef.current);

        renderer.setSetting('nodeReducer', (node, data) => {
            const res = { ...data };
            // colour by the highlighted category the node actually carries, not by
            // its majority vote: that is what makes a rare substrate visible even
            // where it loses the vote in every cluster it appears in
            res.color = paintRef.current(node, data);
            const activeHighlights = highlightRef.current;
            if (activeHighlights.length) {
                const match = bestMatch(data.composition[colorByRef.current] || {}, activeHighlights);
                res.zIndex = match === null ? 0 : 1;
            }

            if (searchHighlightRef.current === node) {
                res.color = searchMarkRef.current;
                res.size = (data.size || 4) * 1.8;
                res.zIndex = 2;
            }
            return res;
        });

        renderer.on('enterNode', ({ node, event }) => {
            setHoverInfo({ x: event.x, y: event.y, attrs: graph.getNodeAttributes(node) });
        });
        renderer.on('leaveNode', () => setHoverInfo(null));
        renderer.on('clickNode', ({ node }) => {
            const attrs = graph.getNodeAttributes(node);
            setHoverInfo(null);
            if (attrs.isCluster) {
                expandClusterRef.current(attrs.clusterId, node);
            }
        });

        return () => {
            resizeObserver.disconnect();
            renderer.kill();
            sigmaRef.current = null;
        };
        // eslint-disable-next-line react-hooks/exhaustive-deps
    }, []);

    // keep refs in sync with state the reducers/handlers need
    useEffect(() => { colorByRef.current = colorBy; setHighlighted([]); refreshSigma(); }, [colorBy, refreshSigma]);
    useEffect(() => { highlightRef.current = highlighted; refreshSigma(); }, [highlighted, refreshSigma]);

    // the full category lists behind the legend's "add" box
    useEffect(() => {
        let cancelled = false;
        Promise.all(CATEGORY_FIELDS.map((field) => fetch(`/api/network/categories?field=${field}`)
            .then((r) => {
                if (!r.ok) throw new Error(`failed to load ${field} categories`);
                return r.json();
            })
            .then((data) => [field, data.categories || []])))
            .then((pairs) => { if (!cancelled) setCategories(Object.fromEntries(pairs)); })
            .catch((err) => console.error(err));
        return () => { cancelled = true; };
    }, []);

    // seed a legend with the most frequent categories the first time we know them;
    // once the user has added or removed anything, leave their selection alone
    useEffect(() => {
        if (!availableCategories.length) return;
        setLegend((prev) => (prev[colorBy] ? prev : {
            ...prev,
            [colorBy]: availableCategories.slice(0, DEFAULT_LEGEND_SIZE).map((c) => c.label),
        }));
    }, [colorBy, availableCategories]);

    const expandClusterRef = useRef(() => {});
    const expandCluster = useCallback(async (clusterId, clusterNodeKey) => {
        const g = graphRef.current;
        if (!g || !g.hasNode(clusterNodeKey)) return;

        try {
            const res = await fetch(`/api/network/cluster/${clusterId}?threshold=${threshold}`);
            if (!res.ok) throw new Error('failed to expand cluster');
            const data = await res.json();

            const basePos = g.getNodeAttributes(clusterNodeKey);
            g.dropNode(clusterNodeKey);

            const clusterKey = `cluster-${clusterId}`;
            const idToKey = new Map();
            data.nodes.forEach((n, i) => {
                const key = `domain-${n.id}`;
                idToKey.set(n.id, key);
                const angle = (2 * Math.PI * i) / Math.max(data.nodes.length, 1);
                const radius = 0.02 + 0.002 * data.nodes.length;
                if (!g.hasNode(key)) {
                    g.addNode(key, {
                        x: basePos.x + Math.cos(angle) * radius,
                        y: basePos.y + Math.sin(angle) * radius,
                        size: 3,
                        label: n.name,
                        isCluster: false,
                        domainId: n.id,
                        clusterId,
                        substrates: n.substrates,
                        genus: n.genus,
                        kingdom: n.kingdom,
                        extendedSignature: n.extended_signature,
                        dominant: {
                            substrate: n.dominant_substrate,
                            genus: n.genus || 'unknown',
                            kingdom: n.kingdom || 'unknown',
                            cluster: clusterKey,
                        },
                        composition: {
                            substrate: compositionOf(n.substrates),
                            genus: { [n.genus || 'unknown']: 1 },
                            kingdom: { [n.kingdom || 'unknown']: 1 },
                            cluster: { [clusterKey]: 1 },
                        },
                    });
                }
            });
            data.edges.forEach((e) => {
                const s = idToKey.get(e.source);
                const t = idToKey.get(e.target);
                if (s && t && s !== t && !g.hasEdge(s, t)) {
                    g.addEdge(s, t, { size: 1, color: edgeColorRef.current, distance: e.distance });
                }
            });

            forceAtlas2.assign(g, { iterations: 80, settings: forceAtlas2.inferSettings(g) });
            refreshSigma();
        } catch (err) {
            console.error(err);
        }
    }, [threshold, refreshSigma]);
    expandClusterRef.current = expandCluster;

    // load (or reload) the top-level cluster graph for a given threshold
    const loadTopLevelGraph = useCallback((t, { isCancelled } = {}) => {
        setLoading(true);
        return fetch(`/api/network/graph?threshold=${t}`)
            .then((r) => {
                if (!r.ok) throw new Error('failed to load network');
                return r.json();
            })
            .then((data) => {
                if (isCancelled?.() || !graphRef.current) return;
                const g = graphRef.current;
                g.clear();

                data.clusters.forEach((c) => {
                    const clusterKey = `cluster-${c.cluster_id}`;
                    g.addNode(clusterKey, {
                        x: Math.random(),
                        y: Math.random(),
                        size: 3 + Math.sqrt(c.size) * 2.2,
                        label: `${c.representative_name}${c.size > 1 ? ` (+${c.size - 1} more)` : ''}`,
                        isCluster: true,
                        clusterId: c.cluster_id,
                        memberCount: c.size,
                        representativeId: c.representative_id,
                        representativeName: c.representative_name,
                        substrateDiversity: c.substrate_diversity,
                        // majority vote of the cluster's members, one vote per domain
                        dominant: {
                            substrate: c.dominant_substrate,
                            genus: c.dominant_genus,
                            kingdom: c.dominant_kingdom,
                            cluster: clusterKey,
                        },
                        composition: {
                            substrate: c.substrate_counts,
                            genus: c.genus_counts,
                            kingdom: c.kingdom_counts,
                            cluster: { [clusterKey]: c.size },
                        },
                    });
                });
                data.edges.forEach((e) => {
                    const s = `cluster-${e.source}`;
                    const t2 = `cluster-${e.target}`;
                    if (g.hasNode(s) && g.hasNode(t2) && !g.hasEdge(s, t2)) {
                        g.addEdge(s, t2, { size: 1, color: edgeColorRef.current, distance: e.distance });
                    }
                });

                forceAtlas2.assign(g, { iterations: 150, settings: forceAtlas2.inferSettings(g) });
                setMeta({ total_domains: data.total_domains, cluster_count: data.clusters.length });
                searchHighlightRef.current = null;
                // cluster ids are only meaningful within one threshold, so a legend of
                // clusters cannot survive a reload, so drop it and let it re-seed
                setLegend((prev) => ({ ...prev, cluster: undefined }));
                // a prior search may have zoomed the camera in on coordinates that no
                // longer mean anything for this freshly laid-out graph, so reset it
                sigmaRef.current?.getCamera().setState({ x: 0.5, y: 0.5, ratio: 1, angle: 0 });
                refreshSigma();
                setLoading(false);
            })
            .catch((err) => {
                console.error(err);
                if (!isCancelled?.()) setLoading(false);
            });
    }, [refreshSigma]);

    // (re)load the top-level cluster graph whenever the threshold changes
    useEffect(() => {
        let cancelled = false;
        setHighlighted([]);
        setNeighborResult(null);
        loadTopLevelGraph(threshold, { isCancelled: () => cancelled });
        return () => { cancelled = true; };
        // eslint-disable-next-line react-hooks/exhaustive-deps
    }, [threshold]);

    // debounced name search
    useEffect(() => {
        if (!searchQuery.trim() || /^[A-Za-z-]{34}$/.test(searchQuery.trim())) {
            setSearchMatches([]);
            return;
        }
        const handle = setTimeout(() => {
            fetch(`/api/network/search_names?q=${encodeURIComponent(searchQuery.trim())}`)
                .then((r) => r.json())
                .then((data) => setSearchMatches(data.matches || []))
                .catch(() => setSearchMatches([]));
        }, 300);
        return () => clearTimeout(handle);
    }, [searchQuery]);

    const locateNode = useCallback(async (nodeKey, clusterId) => {
        const g = graphRef.current;
        if (!g) return;
        if (!g.hasNode(nodeKey) && clusterId !== null && clusterId !== undefined) {
            await expandCluster(clusterId, `cluster-${clusterId}`);
        }
        if (!g.hasNode(nodeKey) || !sigmaRef.current) return;
        searchHighlightRef.current = nodeKey;
        // camera coordinates are sigma's normalized display space, not the raw
        // graph-space x/y stored on the node: sigma.refresh() first so display
        // data reflects the just-added/laid-out node before we read it
        sigmaRef.current.refresh();
        const displayPos = sigmaRef.current.getNodeDisplayData(nodeKey);
        if (displayPos) {
            sigmaRef.current.getCamera().animate({ x: displayPos.x, y: displayPos.y, ratio: 0.15 }, { duration: 600 });
        }
        refreshSigma();
    }, [expandCluster, refreshSigma]);

    const runNeighborSearch = useCallback(async (domainId, signature) => {
        setSearchLoading(true);
        setSearchMatches([]);
        try {
            const params = domainId != null
                ? `domain_id=${domainId}&threshold=${threshold}`
                : `signature=${signature}&threshold=${threshold}`;
            const res = await fetch(`/api/network/neighbors?${params}&k=12`);
            if (!res.ok) throw new Error('search failed');
            const data = await res.json();
            setNeighborResult(data);
            if (domainId != null) {
                await locateNode(`domain-${domainId}`, data.cluster_id);
            }
        } catch (err) {
            console.error(err);
        } finally {
            setSearchLoading(false);
        }
    }, [threshold, locateNode]);

    const handleSearchSubmit = () => {
        const query = searchQuery.trim();
        if (/^[A-Za-z-]{34}$/.test(query)) {
            runNeighborSearch(null, query.toUpperCase());
        }
    };

    const resetView = () => {
        searchHighlightRef.current = null;
        setNeighborResult(null);
        setHighlighted([]);
        loadTopLevelGraph(threshold);
    };

    const modeNoun = colorBy === 'cluster' ? 'cluster' : colorBy;

    const setLegendLabels = (labels) => setLegend((prev) => ({ ...prev, [colorBy]: labels }));

    const addLegendLabel = (label) => {
        // MORE_OPTIONS is the "+N more" row: getOptionDisabled only stops a real
        // pointer (via pointer-events), so reject it where the value is applied
        if (!label || label === MORE_OPTIONS || legendLabels.includes(label)) return;
        setLegendLabels([...legendLabels, label]);
    };

    const removeLegendLabel = (label) => {
        setLegendLabels(legendLabels.filter((l) => l !== label));
        setHighlighted((prev) => prev.filter((l) => l !== label));
    };

    const toggleHighlight = (label) => {
        setHighlighted((prev) => (prev.includes(label) ? prev.filter((l) => l !== label) : [...prev, label]));
    };

    /**
     * Download the graph as it stands (expanded clusters, colours, labels and
     * legend) as a vector SVG plus a PNG of the same figure, zipped together.
     */
    const exportFigure = useCallback(async () => {
        const g = graphRef.current;
        if (!g || !g.order) return;
        setExporting(true);
        try {
            const expanded = g.filterNodes((node, attrs) => !attrs.isCluster).length;
            const subtitle = [
                `${meta.total_domains} domains`,
                `Hamming threshold ${threshold}/34`,
                `coloured by ${modeNoun}`,
                expanded ? `${expanded} domains shown individually` : 'all clusters collapsed',
                highlighted.length ? `highlighting ${highlighted.join(', ')}` : null,
            ].filter(Boolean).join('; ');

            const svg = graphToSvg({
                graph: g,
                colorOf: paintNode,
                showLabels,
                legend: legendLabels.map((label) => ({
                    label,
                    color: colorForLabel(label) || otherColor,
                    count: countByLabel.get(label),
                })),
                caption: { title: 'Sequence similarity network', subtitle },
                colors: {
                    background: surfaceColor,
                    text: theme.palette.text.primary,
                    textSecondary: theme.palette.text.secondary,
                    edge: edgeColor,
                    border: theme.palette.surface.border,
                },
            });

            const png = await svgToPng(svg);
            const stamp = new Date().toISOString().slice(0, 19).replace(/[:T]/g, '-');
            const zip = await createZip([
                { name: 'network.svg', data: svg },
                { name: 'network.png', data: new Uint8Array(await png.arrayBuffer()) },
            ]);
            downloadBlob(zip, `parasect-network-${stamp}.zip`);
        } catch (err) {
            console.error(err);
        } finally {
            setExporting(false);
        }
    }, [meta.total_domains, threshold, highlighted, paintNode, showLabels, legendLabels,
        colorForLabel, otherColor, countByLabel, theme, surfaceColor, edgeColor, modeNoun]);


    return (
        <Box sx={{ display: 'flex', flexDirection: 'column', height: 'calc(100vh - 64px)' }}>
            <Box sx={{ px: { xs: 2, sm: 4 }, py: 2 }}>
                <Typography variant='h4' gutterBottom>
                    Sequence similarity network
                </Typography>
                <Typography variant='body2' color='textSecondary' gutterBottom>
                    Domains from the reference database, grouped into clusters by Hamming distance between their
                    34-residue extended signatures. Click a cluster to expand it into individual domains. A cluster
                    takes the colour its members vote for by majority, one vote per domain; anything outside the
                    legend stays neutral, so every colour on screen is one you can look up.
                    {meta.total_domains > 0 && ` ${meta.total_domains} domains across ${meta.cluster_count} clusters at the current threshold.`}
                </Typography>
            </Box>
            <Divider />

            <Box sx={{ display: 'flex', flexDirection: { xs: 'column', md: 'row' }, flex: 1, minHeight: 0 }}>
                {/* controls sidebar */}
                <Box
                    ref={sidebarRef}
                    sx={{
                        display: sidebar.open ? 'block' : 'none',
                        width: { xs: '100%', md: sidebar.width },
                        flexShrink: 0,
                        p: 2,
                        overflowY: 'auto',
                    }}
                >
                    <Box sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', mb: 1 }}>
                        <Typography variant='subtitle2'>Controls</Typography>
                        <Tooltip title='Hide the controls'>
                            <IconButton size='small' onClick={() => setSidebar((prev) => ({ ...prev, open: false }))}>
                                <KeyboardDoubleArrowLeftIcon fontSize='small' />
                            </IconButton>
                        </Tooltip>
                    </Box>

                    <Typography variant='subtitle2' gutterBottom>
                        Cluster threshold: {displayThreshold} of 34 residues
                    </Typography>
                    <Slider
                        value={displayThreshold}
                        min={0}
                        max={15}
                        step={1}
                        onChange={(e, v) => setDisplayThreshold(v)}
                        onChangeCommitted={(e, v) => setThreshold(v)}
                        valueLabelDisplay='auto'
                        size='small'
                        sx={{ mb: 2 }}
                    />

                    <FormControl size='small' fullWidth sx={{ mb: 2 }}>
                        <InputLabel id='color-by-label'>Color by</InputLabel>
                        <Select
                            labelId='color-by-label'
                            label='Color by'
                            value={colorBy}
                            onChange={(e) => setColorBy(e.target.value)}
                        >
                            <MenuItem value='substrate'>Substrate</MenuItem>
                            <MenuItem value='genus'>Genus</MenuItem>
                            <MenuItem value='kingdom'>Kingdom</MenuItem>
                            <MenuItem value='cluster'>Cluster</MenuItem>
                        </Select>
                    </FormControl>

                    <FormControlLabel
                        control={<Switch size='small' checked={showLabels} onChange={(e) => setShowLabels(e.target.checked)} />}
                        label={<Typography variant='body2'>Node labels</Typography>}
                        sx={{ mb: 0.5, ml: 0 }}
                    />

                    <Box sx={{ display: 'flex', gap: 1, mb: 2 }}>
                        <Button size='small' variant='outlined' fullWidth onClick={resetView}>
                            Reset view
                        </Button>
                        <Button
                            size='small'
                            variant='contained'
                            fullWidth
                            onClick={exportFigure}
                            disabled={exporting || loading}
                            startIcon={exporting
                                ? <CircularProgress size={14} color='inherit' />
                                : <DownloadIcon fontSize='small' />}
                        >
                            {exporting ? 'Packing...' : 'Export'}
                        </Button>
                    </Box>

                    <Divider sx={{ mb: 2 }} />

                    <Box sx={{ display: 'flex', alignItems: 'baseline', justifyContent: 'space-between' }}>
                        <Typography variant='subtitle2'>Legend</Typography>
                        <Typography variant='caption' color='textSecondary'>
                            {legendLabels.length} of {availableCategories.length}
                        </Typography>
                    </Box>
                    <Typography variant='caption' color='textSecondary' display='block' sx={{ mb: 1 }}>
                        Click an entry to highlight everything containing it; everything else greys out.
                    </Typography>

                    <Autocomplete
                        size='small'
                        options={availableCategories.filter((c) => !legendLabels.includes(c.label))}
                        getOptionLabel={(option) => option.label}
                        isOptionEqualToValue={(option, value) => option.label === value.label}
                        filterOptions={(options, state) => {
                            const matches = filterCategories(options, state);
                            if (matches.length <= OPTION_RENDER_LIMIT) return matches;
                            return [
                                ...matches.slice(0, OPTION_RENDER_LIMIT),
                                { label: MORE_OPTIONS, count: matches.length - OPTION_RENDER_LIMIT },
                            ];
                        }}
                        getOptionDisabled={(option) => option.label === MORE_OPTIONS}
                        value={null}
                        blurOnSelect
                        clearOnBlur
                        onChange={(e, option) => addLegendLabel(option?.label)}
                        noOptionsText={availableCategories.length ? 'Nothing left to add' : 'Loading...'}
                        renderOption={(props, option) => (option.label === MORE_OPTIONS ? (
                            <li {...props} key={MORE_OPTIONS} style={{ opacity: 0.7, fontSize: '0.75rem' }}>
                                +{option.count} more - keep typing to narrow
                            </li>
                        ) : (
                            <li {...props} key={option.label}>
                                <Box component='span' sx={{ flex: 1, overflow: 'hidden', textOverflow: 'ellipsis' }}>
                                    {option.label}
                                </Box>
                                <Box component='span' sx={{ ml: 1, color: 'text.secondary', fontSize: '0.75rem' }}>
                                    {option.count}
                                </Box>
                            </li>
                        ))}
                        renderInput={(params) => (
                            <TextField {...params} placeholder={`Add any ${modeNoun}...`} />
                        )}
                        sx={{ mb: 1 }}
                    />

                    <Box sx={{ display: 'flex', flexWrap: 'wrap', gap: 0.5, mb: 1 }}>
                        {legendLabels.map((label) => {
                            const isHighlighted = highlighted.includes(label);
                            const count = countByLabel.get(label);
                            const swatch = colorForLabel(label) || otherColor;
                            return (
                                <Chip
                                    key={label}
                                    size='small'
                                    variant='outlined'
                                    icon={(
                                        <Box sx={{
                                            width: 11, height: 11, borderRadius: '3px', flexShrink: 0,
                                            backgroundColor: swatch, ml: '6px !important',
                                        }} />
                                    )}
                                    label={count === undefined ? label : `${label} (${count})`}
                                    onClick={() => toggleHighlight(label)}
                                    onDelete={() => removeLegendLabel(label)}
                                    sx={{
                                        // the swatch carries the colour and the text stays in
                                        // ink: a name on a saturated fill is the first thing
                                        // that goes unreadable in dark mode
                                        borderColor: swatch,
                                        color: 'text.primary',
                                        opacity: highlighted.length && !isHighlighted ? 0.45 : 1,
                                        fontWeight: isHighlighted ? 700 : 400,
                                        backgroundColor: isHighlighted ? 'action.selected' : 'transparent',
                                    }}
                                />
                            );
                        })}
                        {!legendLabels.length && (
                            <Typography variant='caption' color='textSecondary'>
                                Legend is empty - add a {modeNoun} above, or restore the most frequent ones.
                            </Typography>
                        )}
                    </Box>

                    <Box sx={{ display: 'flex', gap: 1, mb: 1 }}>
                        <Button
                            size='small'
                            onClick={() => setLegendLabels(
                                availableCategories.slice(0, DEFAULT_LEGEND_SIZE).map((c) => c.label),
                            )}
                            disabled={!availableCategories.length}
                        >
                            Top {DEFAULT_LEGEND_SIZE}
                        </Button>
                        <Button
                            size='small'
                            onClick={() => { setLegendLabels([]); setHighlighted([]); }}
                            disabled={!legendLabels.length}
                        >
                            Clear
                        </Button>
                        {highlighted.length > 0 && (
                            <Button size='small' onClick={() => setHighlighted([])}>
                                Unhighlight
                            </Button>
                        )}
                    </Box>

                    {legendLabels.length > SAFE_SLOTS && (
                        <Alert severity='info' variant='outlined' sx={{ py: 0, mb: 1, fontSize: '0.75rem' }}>
                            The first {SAFE_SLOTS} entries use hues checked for red-green colour blindness. Past
                            that they repeat at a second lightness and get harder to tell apart - turn on node
                            labels, or highlight a few at a time.
                        </Alert>
                    )}

                    {matchStats && (
                        <Typography variant='caption' color='textSecondary' display='block' sx={{ mb: 2 }}>
                            {matchStats.nodes === 0
                                ? 'Nothing in the current view contains the highlighted '
                                  + (highlighted.length > 1 ? 'categories' : 'category') + '.'
                                : `${matchStats.nodes} of ${matchStats.total} nodes match `
                                  + `(${matchStats.domains} domain${matchStats.domains === 1 ? '' : 's'}).`}
                        </Typography>
                    )}

                    <Divider sx={{ mb: 2 }} />

                    <Typography variant='subtitle2' gutterBottom>Search</Typography>
                    <TextField
                        size='small'
                        fullWidth
                        placeholder='Domain/protein name or a 34-residue signature'
                        value={searchQuery}
                        onChange={(e) => setSearchQuery(e.target.value)}
                        onKeyDown={(e) => { if (e.key === 'Enter') handleSearchSubmit(); }}
                        sx={{ mb: 1 }}
                    />
                    {searchMatches.length > 0 && (
                        <List dense sx={{ maxHeight: 180, overflowY: 'auto', border: 1, borderColor: 'divider', borderRadius: 1, mb: 1 }}>
                            {searchMatches.map((m) => (
                                <ListItemButton key={m.id} onClick={() => runNeighborSearch(m.id)}>
                                    <ListItemText
                                        primary={m.name}
                                        primaryTypographyProps={{ variant: 'body2', sx: WRAP_ANYWHERE }}
                                    />
                                </ListItemButton>
                            ))}
                        </List>
                    )}
                    {searchLoading && <CircularProgress size={20} />}

                    {neighborResult && (
                        <Paper variant='outlined' sx={{ p: 1.5, mt: 1 }}>
                            <Typography variant='body2' sx={{ fontWeight: 700, ...WRAP_ANYWHERE }}>
                                Nearest neighbors of {neighborResult.query.name || neighborResult.query.signature}
                            </Typography>
                            <Box sx={{ mt: 1, mb: 0.5 }}>
                                <DistanceLegend />
                            </Box>
                            <List dense disablePadding>
                                {neighborResult.neighbors.map((n) => (
                                    <ListItemButton
                                        key={n.id}
                                        onClick={() => runNeighborSearch(n.id)}
                                        sx={{ alignItems: 'flex-start', gap: 1, px: 1, borderRadius: 1 }}
                                    >
                                        <Box sx={{ pt: '2px' }}>
                                            <DistanceBadge distance={n.distance} />
                                        </Box>
                                        <ListItemText
                                            sx={{ my: 0 }}
                                            primary={n.name}
                                            secondary={`${n.substrates.join(', ') || 'unknown substrate'} - ${n.genus}`}
                                            primaryTypographyProps={{ variant: 'body2', sx: WRAP_ANYWHERE }}
                                            secondaryTypographyProps={{ variant: 'caption', sx: WRAP_ANYWHERE }}
                                        />
                                    </ListItemButton>
                                ))}
                            </List>
                        </Paper>
                    )}
                </Box>

                {/* draggable divider; the sigma container's ResizeObserver
                    picks the new canvas size up on its own */}
                {sidebar.open && canResizeSidebar && (
                    <Box
                        role='separator'
                        aria-orientation='vertical'
                        aria-label='Resize the controls'
                        aria-valuenow={sidebar.width}
                        aria-valuemin={SIDEBAR_MIN}
                        aria-valuemax={SIDEBAR_MAX}
                        tabIndex={0}
                        onPointerDown={startSidebarResize}
                        onDoubleClick={() => setSidebarWidth(SIDEBAR_DEFAULT)}
                        onKeyDown={onSidebarResizeKey}
                        sx={{
                            flexShrink: 0,
                            width: '7px',
                            cursor: 'col-resize',
                            // the border is the page's divider line; the handle
                            // is the invisible grab room around it
                            borderLeft: 1,
                            borderColor: 'divider',
                            backgroundColor: 'transparent',
                            transition: 'background-color 120ms ease',
                            '&:hover, &:focus-visible': { backgroundColor: 'primary.main', opacity: 0.35 },
                        }}
                    />
                )}

                {/* graph canvas */}
                <Box sx={{
                    position: 'relative',
                    flex: 1,
                    minWidth: 0,
                    minHeight: 400,
                    borderTop: { xs: 1, md: 0 },
                    borderLeft: { md: sidebar.open && canResizeSidebar ? 0 : 1 },
                    borderColor: 'divider',
                }}>
                    <Box ref={containerRef} sx={{ position: 'absolute', inset: 0 }} />

                    {!sidebar.open && (
                        <Tooltip title='Show the controls'>
                            <IconButton
                                size='small'
                                onClick={() => setSidebar((prev) => ({ ...prev, open: true }))}
                                sx={{
                                    position: 'absolute',
                                    top: 8,
                                    left: 8,
                                    zIndex: 10,
                                    border: 1,
                                    borderColor: 'divider',
                                    backgroundColor: 'background.paper',
                                    '&:hover': { backgroundColor: 'background.paper' },
                                }}
                            >
                                <TuneIcon fontSize='small' />
                            </IconButton>
                        </Tooltip>
                    )}

                    {loading && (
                        <Box sx={{
                            position: 'absolute', inset: 0, display: 'flex', flexDirection: 'column',
                            alignItems: 'center', justifyContent: 'center', backgroundColor: 'background.default', opacity: 0.92,
                        }}>
                            <Loading frame1='paras_loading_1.png' frame2='paras_loading_2.png' />
                        </Box>
                    )}

                    {hoverInfo && (
                        <Paper
                            elevation={3}
                            sx={{
                                position: 'absolute',
                                left: Math.min(hoverInfo.x + 12, (containerRef.current?.clientWidth || 400) - 300),
                                // flip above the cursor near the bottom edge, so the
                                // tooltip never runs off the canvas
                                top: Math.min(hoverInfo.y + 12, (containerRef.current?.clientHeight || 400) - 150),
                                p: 1.5,
                                maxWidth: 300,
                                pointerEvents: 'none',
                                zIndex: 10,
                            }}
                        >
                            <Typography variant='body2' sx={{ fontWeight: 700, ...WRAP_ANYWHERE }}>
                                {hoverInfo.attrs.isCluster && hoverInfo.attrs.memberCount > 1
                                    ? `${hoverInfo.attrs.representativeName} +${hoverInfo.attrs.memberCount - 1} more`
                                    : (hoverInfo.attrs.representativeName || hoverInfo.attrs.label)}
                            </Typography>
                            {hoverInfo.attrs.isCluster ? (
                                <>
                                    <Typography variant='caption' display='block'>
                                        {hoverInfo.attrs.memberCount === 1
                                            ? 'a cluster of one - click to expand'
                                            : `${hoverInfo.attrs.memberCount} domains in this cluster - click to expand`}
                                    </Typography>
                                    <Typography variant='caption' display='block' sx={{ mt: 0.5 }}>
                                        Substrates: {describeComposition(hoverInfo.attrs.composition.substrate)}
                                    </Typography>
                                    <Typography variant='caption' display='block'>
                                        Majority vote: {hoverInfo.attrs.dominant.substrate}
                                    </Typography>
                                    <Typography variant='caption' display='block'>
                                        Genus: {describeComposition(hoverInfo.attrs.composition.genus, 3)}
                                    </Typography>
                                    <Typography variant='caption' display='block'>
                                        Kingdom: {describeComposition(hoverInfo.attrs.composition.kingdom, 3)}
                                    </Typography>
                                </>
                            ) : (
                                <>
                                    <Typography variant='caption' display='block' sx={{ mt: 0.5 }}>
                                        Substrate{(hoverInfo.attrs.substrates || []).length > 1 ? 's' : ''}:{' '}
                                        {(hoverInfo.attrs.substrates || []).join(', ') || 'unknown'}
                                    </Typography>
                                    <Typography variant='caption' display='block'>
                                        Genus: {hoverInfo.attrs.genus}
                                    </Typography>
                                    <Typography variant='caption' display='block'>
                                        Kingdom: {hoverInfo.attrs.kingdom}
                                    </Typography>
                                </>
                            )}
                        </Paper>
                    )}
                </Box>
            </Box>
        </Box>
    );
};

export default NetworkGraph;
