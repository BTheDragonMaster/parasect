import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import { Link as RouterLink } from 'react-router-dom';
import { toast } from 'react-toastify';
import {
  Box,
  Button,
  Card,
  CardActionArea,
  Divider,
  LinearProgress,
  Paper,
  Stack,
  TextField,
  Typography,
  Select,
  MenuItem,
  InputLabel,
} from '@mui/material';
import SearchIcon from '@mui/icons-material/Search';
import ClearIcon from '@mui/icons-material/Clear';
import DownloadIcon from '@mui/icons-material/Download';
import HubIcon from '@mui/icons-material/Hub';
import ArrowForwardIcon from '@mui/icons-material/ArrowForward';
import { DataGrid, GridToolbarContainer, GridPagination } from '@mui/x-data-grid';
import Statistics from '../components/Statistics';
import LazyMultiSelect from '../components/LazyMultiSelect';
import { downloadFile, makeDelimited } from '../utils/tabular';

const DEFAULT_PAGE_SIZE = 100;
const MAX_EXPORT_ROWS = 100000;

// Every preset's Editor reports a plain { paramName: value } object via setParams and
// never SQL text. The server (routes/sql.py, PRESETS) holds the one fixed SQL template
// per preset and binds these as parameters, so user input can never become SQL syntax
// regardless of what's typed (quotes, semicolons, whatever). See /api/sql/preset for
// presets.
const QUERYOPTIONS = [
  {
    key: 'substrate',
    label: 'A-domains by substrate',
    kind: 'preset',
    Editor: ({ params, setParams }) => (
      <LazyMultiSelect
        field='substrate'
        label='Substrates'
        placeholder='Start typing a substrate name...'
        value={params.substrate_name || []}
        onChange={(next) => setParams(next.length ? { substrate_name: next } : {})}
      />
    ),
  },
  {
    key: 'proteinId',
    label: 'Substrate specificities by protein ID',
    kind: 'preset',
    Editor: ({ params, setParams }) => (
      <LazyMultiSelect
        field='protein'
        label='Protein IDs'
        placeholder='e.g. P48633.1'
        value={params.protein_id || []}
        onChange={(next) => setParams(next.length ? { protein_id: next } : {})}
      />
    ),
  },
  {
    key: 'species',
    label: 'Substrate specificities by species',
    kind: 'preset',
    Editor: ({ params, setParams }) => (
      <LazyMultiSelect
        field='species'
        label='Species'
        placeholder='e.g. Streptomyces coelicolor'
        value={params.species || []}
        onChange={(next) => setParams(next.length ? { species: next } : {})}
      />
    ),
  },
  {
    key: 'signature',
    label: 'Substrate specificities by A-domain signature (Hamming <= N)',
    kind: 'preset',
    Editor: ({ params, setParams }) => (
      <Stack direction={{ xs: 'column', sm: 'row' }} spacing={1}>
        <TextField
          label="Signature (max 10)"
          value={params.signature || ''}
          inputProps={{ maxLength: 10 }}
          onChange={(e) => setParams({ ...params, signature: e.target.value })}
          fullWidth
          placeholder="e.g., S/T-A-V-I-G-H-D-L"
        />
        <TextField
          label="Max Hamming distance"
          type="number"
          value={params.max_distance ?? 3}
          inputProps={{ min: 0, max: 10, step: 1 }}
          onChange={(e) => setParams({ ...params, max_distance: e.target.value })}
          sx={{ width: 200 }}
        />
      </Stack>
    ),
  },
];

/** True once a preset's params object has at least one non-empty required value. */
function hasUsableParams(params) {
  return Object.values(params || {}).some((v) => {
    // a multi-select contributes an array, and an empty one is still no filter
    if (Array.isArray(v)) return v.length > 0;
    return v !== '' && v !== null && v !== undefined;
  });
}

const CustomTopToolbar = () => (
  <GridToolbarContainer sx={{ justifyContent: 'flex-end', mb: 1 }}>
    <GridPagination />
  </GridToolbarContainer>
);

const QueryDatabase = () => {
  // which query mode is active
  const [selectedKey, setSelectedKey] = useState(QUERYOPTIONS[0].key);
  const selectedOption = useMemo(
    () => QUERYOPTIONS.find((o) => o.key === selectedKey) || QUERYOPTIONS[0],
    [selectedKey]
  );

  const [presetParams, setPresetParams] = useState({});

  // grid state
  const [rows, setRows] = useState([]);
  const [columns, setColumns] = useState([]);
  const [loading, setLoading] = useState(false);
  const [rowCount, setRowCount] = useState(0);
  const [page, setPage] = useState(0);
  const [pageSize, setPageSize] = useState(DEFAULT_PAGE_SIZE);
  const [sortModel, setSortModel] = useState([]);
  // distinguishes "you haven't searched yet" from "that search found nothing";
  // an empty grid can't say which, so it isn't shown for either
  const [hasSearched, setHasSearched] = useState(false);
  const lastRequestRef = useRef(0);

  // Memoized active sort parameters
  const activeSort = useMemo(() => {
    if (!sortModel.length) return { sortBy: null, sortDir: null };
    const { field, sort } = sortModel[0];
    return { sortBy: field, sortDir: sort };
  }, [sortModel]);

  // Function to build columns from sample rows
  const buildColumnsFromRows = (sampleRows) => {
    if (!sampleRows.length) return [];
    const keys = Object.keys(sampleRows[0]);
    return keys.map((k) => ({
      field: k,
      headerName: k,
      flex: 1,
      minWidth: 120,
    }))
  };

  // Ensure each row has a unique 'id' field for DataGrid
  const ensureRowIds = (arr) => arr.map((r, i) => (r.id ? r : { id: `${page}-${i}`, ...r }));

  // Fetch results from the server via the safe, parameterized preset endpoint
  // (see routes/sql.py: PRESETS). Params are always bound SQL parameters, never
  // interpolated into query text, so this can't be used to inject SQL
  const fetchResults = useCallback(async ({ presetKey, params, p, ps, sortBy, sortDir }) => {
    if (!hasUsableParams(params)) {
      toast.warn('Please fill in the filter above.');
      return;
    }
    setLoading(true);
    const reqId = Date.now();
    lastRequestRef.current = reqId;
    try {
      const res = await fetch('/api/sql/preset', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ preset: presetKey, params, page: p, pageSize: ps, sortBy, sortDir }),
      })
      if (!res.ok) throw new Error(await res.text());
      const data = await res.json();
      if (lastRequestRef.current !== reqId) return;
      const cols = data.columns?.length ? data.columns.map((c) => ({ flex: 1, minWidth: 120, ...c })) : buildColumnsFromRows(data.rows || []);
      const withIds = ensureRowIds(data.rows || []);
      setColumns(cols);
      setRows(withIds);
      setRowCount(Number.isFinite(data.total) ? data.total : withIds.length);
      setHasSearched(true);
    } catch (err) {
      toast.error(`Query failed: ${err.message}`);
      setColumns([]);
      setRows([]);
      setRowCount(0);
      setHasSearched(true);
    } finally {
      setLoading(false);
    }
  }, [page, selectedOption]);

  // Fetch results when page, pageSize, or activeSort changes
  const onSearch = useCallback(() => {
    setPage(0);
    fetchResults({ presetKey: selectedKey, params: presetParams, p: 0, ps: pageSize, ...activeSort });
  }, [selectedKey, presetParams, pageSize, activeSort, fetchResults]);

  // Clear all results and reset state
  const clearAll = () => {
    setRows([]);
    setColumns([]);
    setRowCount(0);
    setPage(0);
    setSortModel([]);
    setHasSearched(false);
    // keep current filter values; if preset is active and its editor cleared, params may be {}
  };

  const exportCurrentPage = (format) => {
    if (!columns.length || !rows.length) {
      toast.info('Nothing to export.');
      return;
    }
    const delim = format === 'csv' ? ',' : '\t';
    const text = makeDelimited(columns, rows, delim);
    downloadFile(text, `results_page${page + 1}.${format}`, format === 'csv' ? 'text/csv' : 'text/tab-separated-values');
  };

  const exportAll = async (format) => {
    setLoading(true);
    try {
      if (!hasUsableParams(presetParams)) return toast.warn('Please fill in the filter above first!');
      const res = await fetch('/api/sql/preset', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          preset: selectedKey,
          params: presetParams,
          page: 0,
          pageSize: MAX_EXPORT_ROWS,
          sortBy: activeSort.sortBy,
          sortDir: activeSort.sortDir,
        }),
      });
      if (!res.ok) throw new Error(await res.text());
      const data = await res.json();
      const cols = data.columns?.length ? data.columns.map((c) => ({ flex: 1, minWidth: 120, ...c })) : buildColumnsFromRows(data.rows || []);
      const delim = format === 'csv' ? ',' : '\t';
      const text = makeDelimited(cols, data.rows || [], delim);
      downloadFile(text, `results_all.${format}`, format === 'csv' ? 'text/csv' : 'text/tab-separated-values');
    } catch (err) {
      toast.error(`Export failed: ${err.message}`);
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    // only refetch if we've already run at least once
    if (!hasUsableParams(presetParams)) return;
    if (rows.length === 0 && rowCount === 0) return;
    fetchResults({ presetKey: selectedKey, params: presetParams, p: page, ps: pageSize, ...activeSort });
  }, [page, pageSize, activeSort.sortBy, activeSort.sortDir]); // eslint-disable-line react-hooks/exhaustive-deps

  return (
    <Box p={2}>
      <Typography
        variant="h5"
        gutterBottom
      >
        Query database
      </Typography>

      <Statistics />

      <Card variant="outlined" sx={{ mt: 2, borderRadius: 3 }}>
        <CardActionArea
          component={RouterLink}
          to="/network"
          sx={{ p: 2, display: 'flex', alignItems: 'center', gap: 2 }}
        >
          <HubIcon color="primary" sx={{ fontSize: 28, flexShrink: 0 }} />
          <Box sx={{ flex: 1, minWidth: 0 }}>
            <Typography variant="subtitle1" sx={{ fontWeight: 700 }}>
              Explore the sequence similarity network
            </Typography>
            <Typography variant="body2" color="text.secondary">
              See how these domains relate to each other, clustered by their signatures, instead of as a table.
            </Typography>
          </Box>
          <ArrowForwardIcon sx={{ flexShrink: 0, display: { xs: 'none', sm: 'block' } }} />
        </CardActionArea>
      </Card>

      <Divider sx={{ my: 2 }} />

      <Stack spacing={1} sx={{ mb: 1 }}>
        <Box sx={{ minWidth: 260 }}>
          <InputLabel id="query-mode-label">Query mode</InputLabel>
          <Select
            fullWidth
            labelId='qtype-label'
            value={selectedKey}
            onChange={(e) => {
              const nextKey = e.target.value;
              setSelectedKey(nextKey);
              // clear filter values and results on mode switch
              setPresetParams({});
              setRows([]);
              setColumns([]);
              setRowCount(0);
              setPage(0);
              setSortModel([]);
              setHasSearched(false);
            }}
          >
            {QUERYOPTIONS.map((o) => (
              <MenuItem key={o.key} value={o.key}>{o.label}</MenuItem>
            ))}
          </Select>
        </Box>

        {selectedOption.Editor ? (
          <selectedOption.Editor params={presetParams} setParams={setPresetParams} />
        ) : null}
      </Stack>

      <Stack
        direction={{ xs: 'column', sm: 'row' }}
        spacing={1}
        alignItems="stretch"
      >
        <Stack spacing={1} minWidth={{ xs: '100%', sm: 220 }}>
          <Stack direction="row" spacing={1}>
            <Button startIcon={<SearchIcon />} variant="contained" onClick={onSearch} disabled={loading} fullWidth>
              Search
            </Button>
            <Button startIcon={<ClearIcon />} variant="outlined" onClick={clearAll} disabled={loading && !rows.length} fullWidth>
              Clear
            </Button>
          </Stack>
          <Stack direction="row" spacing={1}>
            {/* <Button startIcon={<DownloadIcon />} size="small" onClick={() => exportAll('csv')} disabled={loading}>CSV all</Button> */}
            <Button startIcon= {<DownloadIcon />} size="small" onClick={() => exportAll('tsv')} disabled={loading}>Download results</Button>
          </Stack>
        </Stack>
      </Stack>

      <Divider sx={{ my: 2 }} />

      {rows.length > 0 ? (
        <Box sx={{ height: 520, width: '100%' }}>
          <DataGrid
            rows={rows}
            columns={columns}
            pagination
            paginationMode="server"
            sortingMode="server"
            pageSizeOptions={[10, 25, 50, 100]}
            rowCount={rowCount}
            page={page}
            onPaginationModelChange={(model) => {
              if (model.pageSize !== pageSize) setPageSize(model.pageSize);
              if (model.page !== page) setPage(model.page);
            }}
            sortingOrder={["asc", "desc"]}
            sortModel={sortModel}
            onSortModelChange={(model) => setSortModel(model)}
            disableRowSelectionOnClick
            loading={loading}
            slots={{
              toolbar: CustomTopToolbar,
              loadingOverlay: LinearProgress,
            }}
          />
        </Box>
      ) : (
        /* an empty grid is all chrome and no information: a header row, paging
           controls and 500px of nothing. Say what's going on instead. */
        <Paper
          variant="outlined"
          sx={{ py: 5, px: 3, textAlign: 'center', borderStyle: 'dashed' }}
        >
          {loading ? (
            <Typography variant="body2" color="text.secondary">Running your query...</Typography>
          ) : hasSearched ? (
            <>
              <Typography variant="subtitle1" sx={{ fontWeight: 600 }}>
                No matches
              </Typography>
              <Typography variant="body2" color="text.secondary" sx={{ mt: 0.5 }}>
                Nothing in the reference database fits that filter. Try adding more values, or a different query mode.
              </Typography>
            </>
          ) : (
            <>
              <Typography variant="subtitle1" sx={{ fontWeight: 600 }}>
                Pick a filter and hit search
              </Typography>
              <Typography variant="body2" color="text.secondary" sx={{ mt: 0.5 }}>
                You can select several values at once, and every match is combined into one table.
              </Typography>
            </>
          )}
        </Paper>
      )}

      {loading && (
        <Box sx={{ position: 'fixed', left: 0, right: 0, top: 0, zIndex: 1200 }}>
          <LinearProgress />
        </Box>
      )}
    </Box>
  );
};

export default QueryDatabase;
