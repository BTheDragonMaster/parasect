import React, { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import { useSearchParams } from 'react-router-dom';
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
import ViewStreamIcon from '@mui/icons-material/ViewStream';

import Loading from '../components/Loading';
import { DistanceBadge, DistanceLegend } from '../components/DistanceBadge';
import SignatureQueryPanel from '../components/SignatureQueryPanel';
import { useColorMode } from '../theme/ColorModeContext';
import { categoricalColor, MUTED_MARK, OTHER, QUERY_MARK, SAFE_SLOTS } from '../theme';
import { graphToSvg, svgToPng } from '../utils/networkExport';
import {
    NEIGHBOR_COLUMNS, QUERY_COLUMNS, exportStamp, fileSafe, neighborRows, queryRows,
} from '../utils/neighborExport';
import { compareUrl } from '../utils/compare';
import { MAX_QUERIES, checkSignature, describePlacement } from '../utils/signatures';
import { TSV_MIME, downloadFile, makeDelimited } from '../utils/tabular';
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

/** A placed signature expands the clusters it joins, up to this size. Past it
 * (at threshold 10 one cluster holds ~1,500 domains) the signature is linked to
 * the collapsed cluster instead, and clicking the cluster still expands it. */
const AUTO_EXPAND_LIMIT = 300;
const NEIGHBOR_COUNTS = [10, 25, 50];
const QUERY_NODE_SIZE = 7;

/**
 * Signatures handed over in the URL, as the results page's "Show in network"
 * does: repeated signature parameters, each with an optional name at the
 * same position.
 *
 * @param {URLSearchParams} params - the page's query string.
 * @returns {{queries: Array<{name: string|null, signature: string}>, errors: string[]}} -
 *     usable signatures, and one message per unusable one.
 */
function readUrlQueries(params) {
    const names = params.getAll('name');
    const queries = [];
    const errors = [];
    params.getAll('signature').slice(0, MAX_QUERIES).forEach((raw, i) => {
        const name = (names[i] || '').trim() || null;
        const { signature, problem } = checkSignature(raw);
        if (problem) errors.push(`${name || `Signature ${i + 1}`} from the link: ${problem}`);
        else queries.push({ name, signature });
    });
    return { queries, errors };
}

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
    const edgePaintRef = useRef(() => null);
    const themeRef = useRef(null);
    const selectedQueryRef = useRef(null);

    // placed signatures, mirrored for the graph code that runs outside a render
    // (sigma handlers, cluster expansion); state drives the sidebar
    const queriesRef = useRef([]);
    const placementsRef = useRef({});
    const thresholdRef = useRef(5);
    const neighborKRef = useRef(NEIGHBOR_COUNTS[0]);
    // only the latest placement request may touch the graph
    const placeSeqRef = useRef(0);
    const queryCounterRef = useRef(0);
    // settles once the top-level graph for the current threshold is in place
    const graphReadyRef = useRef(Promise.resolve());
    const selectQueryRef = useRef(() => {});
    const resultPanelRef = useRef(null);

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
    const [neighborK, setNeighborK] = useState(NEIGHBOR_COUNTS[0]);

    const [searchParams] = useSearchParams();
    // signatures placed in the graph, in the order they were added
    const [queries, setQueries] = useState([]);
    // query key -> where the server placed it, at the current threshold
    const [placements, setPlacements] = useState({});
    const [selectedQuery, setSelectedQuery] = useState(null);
    const [placing, setPlacing] = useState(false);
    // new signatures only join queriesRef once placed, so a threshold change
    // while they are in flight would re-place the old set and drop them
    const [adding, setAdding] = useState(false);
    const [placeError, setPlaceError] = useState(null);
    const [downloadingResults, setDownloadingResults] = useState(false);

    const [sidebar, setSidebar] = useState(readSidebarPreference);

    thresholdRef.current = threshold;

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
    const queryColor = QUERY_MARK[colorMode];
    // a link to a signature's nearest domain that is beyond the threshold, so
    // it is drawn, but quieter than one that actually joins a cluster
    const queryFaintColor = theme.palette.text.secondary;

    /**
     * The colour a node takes, shared by the renderer and the SVG export so the
     * downloaded figure is the picture on screen.
     */
    const paintNode = useCallback((node, attrs) => {
        // a placed signature is never a category: it keeps its mark through any
        // colour mode or highlight
        if (attrs.isQuery) return queryColor;
        const mode = colorByRef.current;
        const activeHighlights = highlightRef.current;
        if (activeHighlights.length) {
            const match = bestMatch(attrs.composition[mode] || {}, activeHighlights);
            return match === null ? mutedColor : (colorForLabel(match) || otherColor);
        }
        return colorForLabel(attrs.dominant[mode] || 'unknown') || otherColor;
    }, [colorForLabel, mutedColor, otherColor, queryColor]);

    useEffect(() => { paintRef.current = paintNode; sigmaRef.current?.refresh(); }, [paintNode]);

    /** Stroke for a placed signature's links, or null for an ordinary edge; shared like paintNode. */
    const paintEdge = useCallback((edge, attrs) => {
        if (!attrs.queryLink) return null;
        return attrs.beyondThreshold
            ? { color: queryFaintColor, width: 1 }
            : { color: queryColor, width: 1.5 };
    }, [queryColor, queryFaintColor]);

    useEffect(() => { edgePaintRef.current = paintEdge; sigmaRef.current?.refresh(); }, [paintEdge]);

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
        let referenceNodes = 0;
        const activeHighlights = highlightRef.current;
        const mode = colorByRef.current;

        g.forEachNode((node, attrs) => {
            // placed signatures belong to no category and aren't reference domains
            if (attrs.isQuery) return;
            referenceNodes += 1;
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
                ? { nodes: matchedNodes, domains: matchedDomains, total: referenceNodes }
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
            if (data.isQuery) {
                // highlighted keeps the hover ring drawn permanently, on top of
                // everything else; the label is forced through by the attribute
                res.color = paintRef.current(node, data);
                res.highlighted = true;
                if (selectedQueryRef.current === node) res.size = data.size * 1.3;
                return res;
            }
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

        // query links carry their own colour, and the theme effect repaints
        // every edge's attribute, so it is applied here rather than stored
        renderer.setSetting('edgeReducer', (edge, data) => {
            const style = edgePaintRef.current(edge, data);
            return style ? { ...data, color: style.color, size: style.width } : data;
        });

        renderer.on('enterNode', ({ node, event }) => {
            setHoverInfo({ x: event.x, y: event.y, attrs: graph.getNodeAttributes(node) });
        });
        renderer.on('leaveNode', () => setHoverInfo(null));
        renderer.on('clickNode', ({ node }) => {
            const attrs = graph.getNodeAttributes(node);
            setHoverInfo(null);
            if (attrs.isQuery) {
                selectQueryRef.current(attrs.queryKey);
            } else if (attrs.isCluster) {
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

    /**
     * Put every placed signature into the graph and link it to what it joins.
     *
     * Runs after anything that changes the graph (a reload, an expanded cluster,
     * a new placement), so it has to work from whatever is there: each link goes
     * to the domain's own node when its cluster is expanded, and to the
     * collapsed cluster node when it isn't. Signatures that are gone are
     * dropped; ones already in the graph keep their position.
     *
     * @returns {boolean} - whether any signature node was newly added, so needs a layout pass.
     */
    const wireQueries = useCallback(() => {
        const g = graphRef.current;
        if (!g) return false;
        const current = queriesRef.current;
        const byKey = placementsRef.current;
        const wanted = new Set(current.filter((q) => byKey[q.key]).map((q) => q.key));

        g.filterNodes((node, attrs) => attrs.isQuery && !wanted.has(node)).forEach((node) => g.dropNode(node));

        // the jitter below is relative to the layout's own scale, which
        // forceAtlas2 leaves anywhere from ~1 to several hundred units wide
        let minX = Infinity;
        let maxX = -Infinity;
        g.forEachNode((node, attrs) => {
            minX = Math.min(minX, attrs.x);
            maxX = Math.max(maxX, attrs.x);
        });
        const jitter = Number.isFinite(minX) ? Math.max(maxX - minX, 1) * 0.02 : 0.02;

        let added = false;
        current.forEach((q) => {
            const placement = byKey[q.key];
            if (!placement) return;

            // links arrive closest first, so the first one per node is the one to keep
            const targets = new Map();
            placement.links.forEach((link) => {
                const domainKey = `domain-${link.id}`;
                const target = g.hasNode(domainKey) ? domainKey : `cluster-${link.cluster_id}`;
                if (g.hasNode(target) && !targets.has(target)) targets.set(target, link);
            });

            const attrs = {
                label: q.name,
                size: QUERY_NODE_SIZE,
                forceLabel: true,
                isQuery: true,
                queryKey: q.key,
                signature: q.signature,
                // empty rather than missing, so code that walks every node's
                // composition doesn't have to know about signatures
                dominant: {},
                composition: { substrate: {}, genus: {}, kingdom: {}, cluster: {} },
            };
            if (g.hasNode(q.key)) {
                g.edges(q.key).forEach((edge) => g.dropEdge(edge));
                g.mergeNodeAttributes(q.key, attrs);
            } else {
                // start in the middle of whatever it links to, nudged off so it
                // never sits exactly on top of a node
                const points = [...targets.keys()].map((key) => g.getNodeAttributes(key));
                const cx = points.length ? points.reduce((sum, p) => sum + p.x, 0) / points.length : 0;
                const cy = points.length ? points.reduce((sum, p) => sum + p.y, 0) / points.length : 0;
                const angle = Math.random() * 2 * Math.PI;
                g.addNode(q.key, { ...attrs, x: cx + Math.cos(angle) * jitter, y: cy + Math.sin(angle) * jitter });
                added = true;
            }
            targets.forEach((link, target) => {
                g.addEdge(q.key, target, {
                    size: 1,
                    color: edgeColorRef.current,
                    distance: link.distance,
                    queryLink: true,
                    beyondThreshold: link.beyond_threshold,
                });
            });
        });
        return added;
    }, []);

    const expandClusterRef = useRef(() => {});
    const expandCluster = useCallback(async (clusterId, clusterNodeKey, { layout = true } = {}) => {
        const g = graphRef.current;
        if (!g || !g.hasNode(clusterNodeKey)) return;

        try {
            const res = await fetch(`/api/network/cluster/${clusterId}?threshold=${threshold}`);
            if (!res.ok) throw new Error('failed to expand cluster');
            const data = await res.json();
            // the threshold moved while this was in flight: the graph has been
            // rebuilt, and a cluster with this id is now a different cluster
            if (thresholdRef.current !== threshold || !g.hasNode(clusterNodeKey)) return;

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

            // dropping the cluster node took any signature's link to it along;
            // this re-links those signatures to the members themselves
            wireQueries();
            if (layout) forceAtlas2.assign(g, { iterations: 80, settings: forceAtlas2.inferSettings(g) });
            refreshSigma();
        } catch (err) {
            console.error(err);
        }
    }, [threshold, refreshSigma, wireQueries]);
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

                // placed signatures survive a reset: they rejoin the collapsed
                // clusters and settle in the same layout pass. (A threshold change
                // clears their placements first, since cluster ids don't carry over.)
                wireQueries();
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
    }, [refreshSigma, wireQueries]);

    /**
     * Ask server where signatures land at a threshold.
     *
     * @param {Array<{key: string, name: string, signature: string}>} list - signatures to place.
     * @param {number} t - threshold.
     * @returns {Promise<Object<string, object>>} - placement per query key.
     */
    const fetchPlacements = useCallback(async (list, t) => {
        const res = await fetch('/api/network/place', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                threshold: t,
                k: neighborKRef.current,
                queries: list.map(({ name, signature }) => ({ name, signature })),
            }),
        });
        const data = await res.json().catch(() => ({}));
        if (!res.ok) throw new Error(data.error || 'The signatures could not be placed.');
        return Object.fromEntries(list.map((q, i) => [q.key, data.queries[i]]));
    }, []);

    /**
     * Put fresh placements into the graph.
     *
     * @param {Object<string, object>} byKey - placement per query key.
     * @param {object} options - what else to do.
     * @param {string[]} options.expandFor - query keys whose clusters to expand, up to AUTO_EXPAND_LIMIT.
     * @param {() => boolean} options.isStale - true once a newer placement has taken over.
     * @returns {Promise<void>}
     */
    const applyPlacements = useCallback(async (byKey, { expandFor = [], isStale = () => false } = {}) => {
        placementsRef.current = byKey;
        setPlacements(byKey);
        const g = graphRef.current;
        if (!g) return;

        const clusterIds = new Set();
        expandFor.forEach((key) => (byKey[key]?.clusters || []).forEach((c) => {
            if (c.size <= AUTO_EXPAND_LIMIT) clusterIds.add(c.cluster_id);
        }));
        for (const clusterId of clusterIds) {
            const clusterKey = `cluster-${clusterId}`;
            if (g.hasNode(clusterKey)) await expandClusterRef.current(clusterId, clusterKey, { layout: false });
            if (isStale()) return;
        }

        const added = wireQueries();
        if (clusterIds.size || added) {
            forceAtlas2.assign(g, { iterations: 80, settings: forceAtlas2.inferSettings(g) });
        }
        refreshSigma();
    }, [wireQueries, refreshSigma]);

    /** Pan to a node, and optionally give it the search marker. */
    const focusNode = useCallback((nodeKey, { mark = true } = {}) => {
        const g = graphRef.current;
        if (!g || !g.hasNode(nodeKey) || !sigmaRef.current) return;
        searchHighlightRef.current = mark ? nodeKey : null;
        // camera coordinates are sigma's normalized display space, not the raw
        // graph-space x/y stored on the node: sigma.refresh() first so display
        // data reflects the just-added/laid-out node before we read it
        sigmaRef.current.refresh();
        const displayPos = sigmaRef.current.getNodeDisplayData(nodeKey);
        if (displayPos) {
            sigmaRef.current.getCamera().animate({ x: displayPos.x, y: displayPos.y, ratio: 0.15 }, { duration: 600 });
        }
        refreshSigma();
    }, [refreshSigma]);

    /** List a placed signature's neighbours in the sidebar, and optionally pan to it. */
    const selectQuery = useCallback((key, { locate = true } = {}) => {
        selectedQueryRef.current = key;
        setSelectedQuery(key);
        setNeighborResult(null);
        if (key && locate) focusNode(key, { mark: false });
        else sigmaRef.current?.refresh();
    }, [focusNode]);
    selectQueryRef.current = selectQuery;

    /**
     * Place a set of signatures, replacing whatever placements there were.
     *
     * @param {Array<{key: string, name: string, signature: string}>} next - every signature to show.
     * @param {object} options - what else to do.
     * @param {string[]} options.expandFor - query keys whose clusters to expand.
     * @param {string|null} options.focusKey - query to select and pan to afterwards.
     * @returns {Promise<boolean>} - whether they were placed.
     */
    const placeAll = useCallback(async (next, { expandFor = [], focusKey = null } = {}) => {
        placeSeqRef.current += 1;
        const seq = placeSeqRef.current;
        setPlacing(true);
        setPlaceError(null);
        try {
            const byKey = next.length ? await fetchPlacements(next, thresholdRef.current) : {};
            await graphReadyRef.current;
            if (seq !== placeSeqRef.current) return false;
            queriesRef.current = next;
            setQueries(next);
            const isStale = () => seq !== placeSeqRef.current;
            await applyPlacements(byKey, { expandFor, isStale });
            if (focusKey && !isStale()) selectQuery(focusKey);
            return true;
        } catch (err) {
            console.error(err);
            if (seq === placeSeqRef.current) setPlaceError(err.message);
            return false;
        } finally {
            if (seq === placeSeqRef.current) setPlacing(false);
        }
    }, [fetchPlacements, applyPlacements, selectQuery]);

    /**
     * Add signatures to the ones already placed.
     *
     * @param {Array<{name: string|null, signature: string}>} parsed - checked signatures.
     * @returns {Promise<boolean>} - whether they are now in the graph.
     */
    const addQueries = useCallback(async (parsed) => {
        const existing = queriesRef.current;
        const fresh = [];
        parsed.forEach((p) => {
            const duplicate = [...existing, ...fresh].some(
                (q) => q.signature === p.signature && (!p.name || q.name === p.name),
            );
            if (duplicate) return;
            queryCounterRef.current += 1;
            const n = queryCounterRef.current;
            fresh.push({ key: `query-${n}`, name: p.name || `signature_${n}`, signature: p.signature });
        });

        if (!fresh.length) {
            // everything pasted is already there: point at it instead
            const match = existing.find((q) => q.signature === parsed[0]?.signature);
            if (match) selectQuery(match.key);
            return true;
        }
        if (existing.length + fresh.length > MAX_QUERIES) {
            setPlaceError(`Up to ${MAX_QUERIES} signatures can be placed at once, and this would make `
                + `${existing.length + fresh.length}. Remove some first.`);
            return false;
        }
        setAdding(true);
        try {
            return await placeAll([...existing, ...fresh], {
                expandFor: fresh.map((q) => q.key),
                focusKey: fresh[0].key,
            });
        } finally {
            setAdding(false);
        }
    }, [placeAll, selectQuery]);

    const removeQuery = useCallback((key) => {
        const next = queriesRef.current.filter((q) => q.key !== key);
        const { [key]: removed, ...rest } = placementsRef.current;
        queriesRef.current = next;
        placementsRef.current = rest;
        setQueries(next);
        setPlacements(rest);
        if (selectedQueryRef.current === key) selectQuery(null);
        wireQueries();
        refreshSigma();
    }, [selectQuery, wireQueries, refreshSigma]);

    const clearQueries = useCallback(() => {
        queriesRef.current = [];
        placementsRef.current = {};
        setQueries([]);
        setPlacements({});
        setPlaceError(null);
        selectQuery(null);
        wireQueries();
        refreshSigma();
    }, [selectQuery, wireQueries, refreshSigma]);

    // (re)load the top-level cluster graph whenever the threshold changes, and
    // re-place any signatures: which clusters they join depends on it
    useEffect(() => {
        let cancelled = false;
        setHighlighted([]);
        setNeighborResult(null);
        // cluster ids from the old threshold would link signatures to the wrong
        // clusters in the new graph, so they come out until re-placed
        placementsRef.current = {};
        setPlacements({});
        const ready = loadTopLevelGraph(threshold, { isCancelled: () => cancelled });
        graphReadyRef.current = ready;
        if (queriesRef.current.length) {
            const current = queriesRef.current;
            placeAll(current, { expandFor: current.map((q) => q.key), focusKey: selectedQueryRef.current });
        }
        return () => { cancelled = true; };
        // eslint-disable-next-line react-hooks/exhaustive-deps
    }, [threshold]);

    // signatures handed over in the URL (the results page's "Show in network");
    // placeAll waits for the first graph load on its own
    useEffect(() => {
        const { queries: fromUrl, errors } = readUrlQueries(searchParams);
        if (fromUrl.length) addQueries(fromUrl).then(() => { if (errors.length) setPlaceError(errors.join(' ')); });
        else if (errors.length) setPlaceError(errors.join(' '));
        // eslint-disable-next-line react-hooks/exhaustive-deps
    }, []);

    // bring the neighbour list into view when a new one opens: the sidebar is
    // usually scrolled to the controls at the top when a node is clicked
    const activeResultKey = selectedQuery || neighborResult?.query?.id || neighborResult?.query?.signature || null;
    useEffect(() => {
        if (activeResultKey) resultPanelRef.current?.scrollIntoView({ block: 'nearest', behavior: 'smooth' });
    }, [activeResultKey]);

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

    /** Pan to a domain, expanding its cluster first if it is collapsed. */
    const locateNode = useCallback(async (nodeKey, clusterId) => {
        const g = graphRef.current;
        if (!g) return;
        if (!g.hasNode(nodeKey) && clusterId !== null && clusterId !== undefined) {
            await expandCluster(clusterId, `cluster-${clusterId}`);
        }
        focusNode(nodeKey);
    }, [expandCluster, focusNode]);

    const runNeighborSearch = useCallback(async (domainId) => {
        setSearchLoading(true);
        setSearchMatches([]);
        try {
            const res = await fetch(
                `/api/network/neighbors?domain_id=${domainId}&threshold=${threshold}&k=${neighborKRef.current}`,
            );
            if (!res.ok) throw new Error('search failed');
            const data = await res.json();
            selectQuery(null, { locate: false });
            setNeighborResult(data);
            await locateNode(`domain-${domainId}`, data.cluster_id);
        } catch (err) {
            console.error(err);
        } finally {
            setSearchLoading(false);
        }
    }, [threshold, locateNode, selectQuery]);

    const handleSearchSubmit = async () => {
        const query = searchQuery.trim();
        // a signature typed into search gets placed like a pasted one, so it
        // shows up in the graph rather than only as a list
        if (/^[A-Za-z-]{34}$/.test(query) && await addQueries([{ name: null, signature: query.toUpperCase() }])) {
            setSearchQuery('');
        }
    };

    const changeNeighborK = (k) => {
        neighborKRef.current = k;
        setNeighborK(k);
        if (neighborResult?.query?.id != null) runNeighborSearch(neighborResult.query.id);
        // links don't depend on k, only the neighbour lists do, so nothing expands
        if (queriesRef.current.length) placeAll(queriesRef.current);
    };

    const resetView = () => {
        searchHighlightRef.current = null;
        setNeighborResult(null);
        setHighlighted([]);
        graphReadyRef.current = loadTopLevelGraph(threshold);
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
     * The graph as it stands (expanded clusters, placed signatures, colours,
     * labels and legend) as a vector SVG plus a PNG of the same figure.
     *
     * @returns {Promise<Array<{name: string, data: (string|Uint8Array)}>>} - zip entries.
     */
    const renderFigureFiles = useCallback(async () => {
        const g = graphRef.current;
        const expanded = g.filterNodes((node, attrs) => !attrs.isCluster && !attrs.isQuery).length;
        const placed = g.filterNodes((node, attrs) => attrs.isQuery).length;
        const subtitle = [
            `${meta.total_domains} domains`,
            `Hamming threshold ${threshold}/34`,
            `coloured by ${modeNoun}`,
            expanded ? `${expanded} domains shown individually` : 'all clusters collapsed',
            placed ? `${placed} placed signature${placed === 1 ? '' : 's'}` : null,
            highlighted.length ? `highlighting ${highlighted.join(', ')}` : null,
        ].filter(Boolean).join('; ');

        const svg = graphToSvg({
            graph: g,
            colorOf: paintNode,
            edgeStyleOf: paintEdge,
            ringOf: (node, attrs) => (attrs.isQuery ? queryColor : null),
            showLabels,
            legend: [
                ...(placed ? [{ label: 'Your signatures', color: queryColor, ring: true }] : []),
                ...legendLabels.map((label) => ({
                    label,
                    color: colorForLabel(label) || otherColor,
                    count: countByLabel.get(label),
                })),
            ],
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
        return [
            { name: 'network.svg', data: svg },
            { name: 'network.png', data: new Uint8Array(await png.arrayBuffer()) },
        ];
    }, [meta.total_domains, threshold, highlighted, paintNode, paintEdge, queryColor, showLabels, legendLabels,
        colorForLabel, otherColor, countByLabel, theme, surfaceColor, edgeColor, modeNoun]);

    /** Download the figure on its own, as a zip of the SVG and PNG. */
    const exportFigure = useCallback(async () => {
        const g = graphRef.current;
        if (!g || !g.order) return;
        setExporting(true);
        try {
            const zip = await createZip(await renderFigureFiles());
            downloadBlob(zip, `parasect-network-${exportStamp()}.zip`);
        } catch (err) {
            console.error(err);
        } finally {
            setExporting(false);
        }
    }, [renderFigureFiles]);

    /**
     * Download everything about the placed signatures: where each landed, all
     * their neighbours, the figure they appear in, and the settings behind it.
     */
    const downloadQueryResults = useCallback(async () => {
        const placed = queries.filter((q) => placements[q.key]);
        if (!placed.length) return;
        setDownloadingResults(true);
        try {
            const neighbors = placed.flatMap((q) => neighborRows({
                name: q.name, signature: q.signature, threshold, neighbors: placements[q.key].neighbors,
            }));
            const parameters = {
                threshold,
                neighbors_per_signature: neighborK,
                reference_domains: meta.total_domains,
                clusters_at_threshold: meta.cluster_count,
                exported: new Date().toISOString(),
                app_version: process.env.REACT_APP_VERSION || null,
                note: 'Distances are Hamming distances between 34-residue extended signatures. '
                    + 'Cluster ids are only meaningful at the threshold above.',
            };
            const zip = await createZip([
                { name: 'signatures.tsv', data: makeDelimited(QUERY_COLUMNS, queryRows(placed, placements, threshold), '\t') },
                { name: 'neighbors.tsv', data: makeDelimited(NEIGHBOR_COLUMNS, neighbors, '\t') },
                ...await renderFigureFiles(),
                { name: 'parameters.json', data: JSON.stringify(parameters, null, 2) },
            ]);
            downloadBlob(zip, `parasect-signatures-${exportStamp()}.zip`);
        } catch (err) {
            console.error(err);
            setPlaceError('The download could not be built.');
        } finally {
            setDownloadingResults(false);
        }
    }, [queries, placements, threshold, neighborK, meta, renderFigureFiles]);

    // the neighbour list on screen: a placed signature's, or a database domain's
    const selectedPlacement = selectedQuery ? placements[selectedQuery] : null;
    const selectedQueryInfo = queries.find((q) => q.key === selectedQuery);
    let activeResult = null;
    if (selectedPlacement && selectedQueryInfo) {
        activeResult = {
            isQuery: true,
            name: selectedQueryInfo.name,
            signature: selectedQueryInfo.signature,
            neighbors: selectedPlacement.neighbors,
            placement: selectedPlacement,
        };
    } else if (neighborResult) {
        activeResult = {
            isQuery: false,
            id: neighborResult.query.id,
            name: neighborResult.query.name || neighborResult.query.signature,
            signature: neighborResult.query.signature,
            neighbors: neighborResult.neighbors,
        };
    }
    
    const openInCompare = () => {
        if (!activeResult) return;
        const neighborIds = activeResult.neighbors.map((n) => n.id);
        const url = activeResult.isQuery || activeResult.id === null || activeResult.id === undefined
            ? compareUrl({ custom: [{ name: activeResult.name, signature: activeResult.signature }], ids: neighborIds })
            : compareUrl({ ids: [activeResult.id, ...neighborIds] });
        window.open(url, '_blank', 'noopener');
    };

    /** Download the neighbour list on screen as a TSV. */
    const downloadNeighborTable = () => {
        if (!activeResult) return;
        const rows = neighborRows({ ...activeResult, threshold });
        downloadFile(
            makeDelimited(NEIGHBOR_COLUMNS, rows, '\t'),
            `parasect-neighbors-${fileSafe(activeResult.name)}-${exportStamp()}.tsv`,
            TSV_MIME,
        );
    };


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
                    legend stays neutral, so every colour on screen is one you can look up. Paste your own extended
                    signatures under "Your signatures" to see which clusters they would join.
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
                        disabled={adding}
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
                                +{option.count} more. Keep typing to narrow!
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

                    <SignatureQueryPanel
                        queries={queries}
                        placements={placements}
                        selectedKey={selectedQuery}
                        threshold={threshold}
                        placing={placing}
                        error={placeError}
                        downloading={downloadingResults}
                        queryColor={queryColor}
                        onPlace={addQueries}
                        onSelect={(key) => selectQuery(key)}
                        onRemove={removeQuery}
                        onClear={clearQueries}
                        onDownload={downloadQueryResults}
                    />

                    <Divider sx={{ my: 2 }} />

                    <Typography variant='subtitle2' gutterBottom>Search the database</Typography>
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

                    {activeResult && (
                        <Paper ref={resultPanelRef} variant='outlined' sx={{ p: 1.5, mt: 1, scrollMarginBottom: 16 }}>
                            <Box sx={{ display: 'flex', alignItems: 'flex-start', gap: 0.5 }}>
                                <Typography variant='body2' sx={{ fontWeight: 700, flex: 1, ...WRAP_ANYWHERE }}>
                                    Nearest neighbors of {activeResult.name}
                                </Typography>
                                <Select
                                    variant='standard'
                                    disableUnderline
                                    value={neighborK}
                                    onChange={(e) => changeNeighborK(e.target.value)}
                                    renderValue={(k) => `top ${k}`}
                                    inputProps={{ 'aria-label': 'Neighbours to list' }}
                                    sx={{ fontSize: '0.8125rem', flexShrink: 0 }}
                                >
                                    {NEIGHBOR_COUNTS.map((k) => <MenuItem key={k} value={k}>top {k}</MenuItem>)}
                                </Select>
                                <Tooltip title='Compare these signatures side by side (new tab)'>
                                    <IconButton size='small' onClick={openInCompare} sx={{ mt: -0.5 }}>
                                        <ViewStreamIcon fontSize='small' />
                                    </IconButton>
                                </Tooltip>
                                <Tooltip title='Download this list (TSV)'>
                                    <IconButton size='small' onClick={downloadNeighborTable} sx={{ mt: -0.5 }}>
                                        <DownloadIcon fontSize='small' />
                                    </IconButton>
                                </Tooltip>
                            </Box>
                            {activeResult.isQuery && (
                                <>
                                    <Typography
                                        variant='caption'
                                        display='block'
                                        sx={{ fontFamily: 'monospace', color: 'text.secondary', ...WRAP_ANYWHERE }}
                                    >
                                        {activeResult.signature}
                                    </Typography>
                                    <Typography variant='caption' display='block' sx={{ mt: 0.5, ...WRAP_ANYWHERE }}>
                                        {describePlacement(activeResult.placement, threshold)}
                                    </Typography>
                                </>
                            )}
                            <Box sx={{ mt: 1, mb: 0.5 }}>
                                <DistanceLegend />
                            </Box>
                            <List dense disablePadding>
                                {activeResult.neighbors.map((n) => (
                                    <ListItemButton
                                        key={n.id}
                                        // a signature's list stays put and the click just finds the
                                        // domain; a database domain's list walks on to the neighbour's
                                        onClick={() => (activeResult.isQuery
                                            ? locateNode(`domain-${n.id}`, n.cluster_id)
                                            : runNeighborSearch(n.id))}
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
                            {hoverInfo.attrs.isQuery && (
                                <>
                                    <Typography variant='caption' display='block'>
                                        Your signature - click to list its neighbours
                                    </Typography>
                                    <Typography
                                        variant='caption'
                                        display='block'
                                        sx={{ fontFamily: 'monospace', mt: 0.5, ...WRAP_ANYWHERE }}
                                    >
                                        {hoverInfo.attrs.signature}
                                    </Typography>
                                    <Typography variant='caption' display='block' sx={{ mt: 0.5 }}>
                                        {describePlacement(placements[hoverInfo.attrs.queryKey], threshold)}
                                    </Typography>
                                </>
                            )}
                            {hoverInfo.attrs.isQuery ? null : hoverInfo.attrs.isCluster ? (
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
