import React, { useEffect, useMemo, useRef, useState } from 'react';
import {
  Box,
  Card,
  CardContent,
  CardHeader,
  LinearProgress,
  Tab,
  Tabs,
  Typography,
  useMediaQuery,
  useTheme,
} from '@mui/material';
import {
  Bar,
  BarChart,
  Cell,
  LabelList,
  Legend,
  Pie,
  PieChart,
  ResponsiveContainer,
  Tooltip,
  XAxis,
  YAxis,
} from 'recharts';
import { toast } from 'react-toastify';

import { useColorMode } from '../theme/ColorModeContext';
import { categoricalColor, OTHER } from '../theme';

/**
 * The breakdowns you can page through, in tab order.
 */
const VIEWS = [
  { key: 'domain', tab: 'Taxonomic domain', noun: 'taxonomic domain', plural: 'taxonomic domains', chart: 'pie' },
  { key: 'phylum', tab: 'Phylum', noun: 'phylum', plural: 'phyla', chart: 'bar' },
  { key: 'genus', tab: 'Genus', noun: 'genus', plural: 'genera', chart: 'bar' },
  { key: 'species', tab: 'Species', noun: 'species', plural: 'species', chart: 'bar' },
  { key: 'substrate', tab: 'Substrate', noun: 'substrate', plural: 'substrates', chart: 'bar' },
];

const SUBHEADER = 'Counts of distinct adenylation domains that have substrate annotations';

/**
 * A pie puts every slice next to every other, so identity has to survive the
 * all-pairs check rather than just neighbouring pairs, and that caps the
 * validated palette at three hues. Anything past the third-largest taxonomic
 * domain is a sliver anyway, so it folds into one "other" slice instead of
 * spending a colour nobody can distinguish.
 */
const MAX_SLICES = 3;

/** How many bars the ranked views show before the caption takes over. */
const TOP_BARS = 12;

/** Plot height, the same for every view so switching tabs doesn't move the page. */
const PLOT_HEIGHT = 340;

const formatCount = (n) => Number(n).toLocaleString();

export default function Statistics() {
  const [viewKey, setViewKey] = useState(VIEWS[0].key);
  // the rows carry the view they belong to, so a slow response for the tab you
  // just left can't be drawn under the heading of the one you're on now
  const [loaded, setLoaded] = useState({ key: null, rows: [] });
  const [loading, setLoading] = useState(false);
  const cacheRef = useRef(new Map());
  const { mode: colorMode } = useColorMode();
  const theme = useTheme();
  const narrow = useMediaQuery(theme.breakpoints.down('sm'));

  const view = useMemo(() => VIEWS.find((v) => v.key === viewKey) || VIEWS[0], [viewKey]);

  // Fetch a view the first time it's opened; tabs you come back to are already here
  useEffect(() => {
    const cached = cacheRef.current.get(viewKey);
    if (cached) {
      setLoaded({ key: viewKey, rows: cached });
      return undefined;
    }
    let cancelled = false;
    setLoading(true);
    (async () => {
      try {
        const res = await fetch(`/api/sql/stats?view=${encodeURIComponent(viewKey)}`);
        if (!res.ok) throw new Error(await res.text());
        const data = await res.json();
        const received = Array.isArray(data.rows) ? data.rows : [];
        cacheRef.current.set(viewKey, received);
        if (!cancelled) setLoaded({ key: viewKey, rows: received });
      } catch (e) {
        // mark the view as settled-but-failed, otherwise the card sits on
        // "Loading..." for good once the toast has come and gone. Failures aren't
        // cached, so coming back to the tab tries again.
        if (!cancelled) {
          setLoaded({ key: viewKey, rows: [], failed: true });
          toast.error(`Failed to load stats: ${e.message}`);
        }
      } finally {
        if (!cancelled) setLoading(false);
      }
    })();
    return () => { cancelled = true; };
  }, [viewKey]);

  // null until this view's own rows are in which is distinct from "loaded, but empty"
  const rows = loaded.key === viewKey ? loaded.rows : null;

  const categories = useMemo(() => (rows ?? [])
    .map((r) => ({ name: r.label ?? 'unclassified', value: Number(r.value ?? 0) }))
    .filter((d) => d.value > 0)
    .sort((a, b) => b.value - a.value), [rows]);

  const pieData = useMemo(() => {
    const head = categories.slice(0, MAX_SLICES);
    const tail = categories.slice(MAX_SLICES);
    const slices = head.map((d, i) => ({ ...d, color: categoricalColor(i, colorMode) }));
    if (tail.length) {
      slices.push({
        name: `other (${tail.length})`,
        value: tail.reduce((sum, d) => sum + d.value, 0),
        color: OTHER[colorMode],
      });
    }
    return slices.map((d) => ({ ...d, label: `${d.name} (${formatCount(d.value)})` }));
  }, [categories, colorMode]);

  // recharts lays a vertical-layout category axis out in data order, top down,
  // so the server's largest-first order is already the ranking a reader expects
  const barData = useMemo(() => categories.slice(0, TOP_BARS), [categories]);

  // How many categories the chart actually draws, so the caption can be honest
  // about the tail rather than leaving the reader to assume there isn't one
  const shown = view.chart === 'pie' ? MAX_SLICES : TOP_BARS;
  const caption = rows === null || categories.length === 0
    ? ''
    : categories.length > shown
      ? `Showing the ${shown} largest of ${formatCount(categories.length)} ${view.plural}`
        + (view.chart === 'pie' ? '; the rest are grouped as "other".' : '.')
      : `All ${formatCount(categories.length)} ${categories.length === 1 ? view.noun : view.plural}.`;

  // recharts measures tick text at the default font size rather than the 12px it
  // draws at, so a long binomial ("Pseudomonas syringae pv. syringae") wraps to a
  // second line well before it actually fills the axis. Two lines is the better
  // outcome anyway, and a clipped species name is worse than a tall row, so the
  // cap is only here to bound how many lines that can become. A small device has no room
  // for either, so there names are clipped and the tooltip carries them in full.
  const axisWidth = narrow ? 116 : 250;
  const maxLabelChars = narrow ? 13 : 44;

  const tooltipStyle = {
    backgroundColor: theme.palette.background.paper,
    border: `1px solid ${theme.palette.divider}`,
    borderRadius: 8,
    color: theme.palette.text.primary,
  };

  // recharts prints the value in the mark's own colour; the swatch beside it
  // already carries that, so the reading itself stays in ink
  const tooltipItemStyle = { color: theme.palette.text.primary };

  return (
    <Card>
      <CardHeader
        title={`Annotated A-domains per ${view.noun}`}
        subheader={SUBHEADER}
      />
      <Tabs
        value={viewKey}
        onChange={(_, next) => setViewKey(next)}
        variant="scrollable"
        scrollButtons="auto"
        allowScrollButtonsMobile
        aria-label="Breakdown to chart"
        sx={{ px: 2, borderBottom: 1, borderColor: 'divider', minHeight: 40 }}
      >
        {VIEWS.map((v) => (
          <Tab key={v.key} value={v.key} label={v.tab} sx={{ minHeight: 40, py: 0 }} />
        ))}
      </Tabs>
      <CardContent sx={{ position: 'relative' }}>
        {loading && (
          <Box sx={{ position: 'absolute', left: 0, right: 0, top: 0 }}>
            <LinearProgress />
          </Box>
        )}
        {rows === null || categories.length === 0 ? (
          <Box sx={{ height: PLOT_HEIGHT, display: 'flex', alignItems: 'center', justifyContent: 'center' }}>
            <Typography variant="body2" color="text.secondary">
              {rows === null ? 'Loading...' : loaded.failed ? 'Could not load this breakdown' : 'No data'}
            </Typography>
          </Box>
        ) : view.chart === 'pie' ? (
          <ResponsiveContainer width="100%" height={PLOT_HEIGHT}>
            <PieChart>
              <Pie
                data={pieData}
                dataKey="value"
                nameKey="label"
                outerRadius="80%"
                isAnimationActive={false}
              >
                {pieData.map((entry) => (
                  <Cell
                    key={entry.name}
                    fill={entry.color}
                    // a 2px gap in the surface colour keeps touching slices apart
                    // for anyone the hues alone don't separate
                    stroke={theme.palette.background.paper}
                    strokeWidth={2}
                  />
                ))}
              </Pie>
              <Tooltip contentStyle={tooltipStyle} itemStyle={tooltipItemStyle} />
              <Legend
                // recharts paints legend text in the slice colour by default; the
                // swatch already carries identity, so the words stay in ink
                formatter={(value) => (
                  <span style={{ color: theme.palette.text.primary }}>{value}</span>
                )}
              />
            </PieChart>
          </ResponsiveContainer>
        ) : (
          <ResponsiveContainer width="100%" height={PLOT_HEIGHT}>
            <BarChart data={barData} layout="vertical" margin={{ top: 4, right: 56, bottom: 4, left: 4 }}>
              {/* one bar per category and the count written at every bar end, so
                  the hue carries no information. Every bar gets slot 1, and an
                  x-axis would only repeat what the labels already say */}
              <XAxis type="number" hide domain={[0, 'dataMax']} />
              <YAxis
                type="category"
                dataKey="name"
                width={axisWidth}
                interval={0}
                tickLine={false}
                axisLine={false}
                tick={{ fill: theme.palette.text.secondary, fontSize: 12 }}
                // trimmed before the ellipsis so the clipped label is one
                // unbroken word -- recharts only wraps on spaces, and
                // "Pseudomonas ..." would otherwise drop the ... onto its own line
                tickFormatter={(v) => {
                  const label = String(v);
                  return label.length > maxLabelChars
                    ? `${label.slice(0, maxLabelChars - 3).trimEnd()}...`
                    : label;
                }}
              />
              <Tooltip
                cursor={{ fill: theme.palette.action.hover }}
                contentStyle={tooltipStyle}
                itemStyle={tooltipItemStyle}
                formatter={(value) => [formatCount(value), 'A-domains']}
              />
              <Bar
                dataKey="value"
                fill={categoricalColor(0, colorMode)}
                barSize={14}
                radius={[0, 4, 4, 0]}
                isAnimationActive={false}
              >
                <LabelList
                  dataKey="value"
                  position="right"
                  formatter={formatCount}
                  fill={theme.palette.text.secondary}
                  fontSize={12}
                />
              </Bar>
            </BarChart>
          </ResponsiveContainer>
        )}
        <Typography
          variant="caption"
          color="text.secondary"
          // reserved whether or not there's a tail to describe, so the card is
          // the same height on every tab
          sx={{ display: 'block', minHeight: 20, mt: 1, textAlign: 'center' }}
        >
          {caption}
        </Typography>
      </CardContent>
    </Card>
  );
}
