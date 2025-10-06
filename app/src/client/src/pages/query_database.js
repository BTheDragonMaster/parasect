import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import { toast } from 'react-toastify';
import { 
  Box, 
  Button,
  Divider,
  LinearProgress,
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
import { DataGrid, GridToolbarContainer, GridPagination } from '@mui/x-data-grid';
import Statistics from '../components/Statistics';

const DEFAULT_PAGE_SIZE = 100;
const MAX_EXPORT_ROWS = 100000;

const QUERYOPTIONS = [
  // {
  //   key: 'free',
  //   label: 'Free-form SQL',
  //   kind: 'free',
  //   defaultQuery: 'SELECT name, smiles FROM substrate LIMIT 500',
  //   Editor: ({ query, setQuery }) => (
  //     <TextField
  //       label="SQL query"
  //       value={query}
  //       onChange={(e) => setQuery(e.target.value)}
  //       fullWidth
  //       multiline
  //       minRows={3}
  //       placeholder="e.g., SELECT * FROM substrate LIMIT 500"
  //     />
  //   ),
  // },
  {
    key: 'proteinId',
    label: 'Substrate specificities by protein ID',
    kind: 'preset',
    defaultQuery: "SELECT * FROM protein_synonym LIMIT 500",
    // Editor compiles SQL directly into the main query state
    Editor: ({ presetInput, setPresetInput, setQuery }) => (
      <TextField
        label="Protein ID"
        value={presetInput}
        onChange={(e) => {
          const v = e.target.value;
          setPresetInput(v);
          const t = v.trim();
          setQuery(
            t
              ? `
    SELECT
      ad.id                    AS domain_id,
      p.id                     AS protein_id,
      ps.synonym               AS protein_synonym,
      pda.domain_number,
      ad.signature,
      ad.extended_signature,
      s.name                   AS substrate_name,
      s.smiles                 AS substrate_smiles
    FROM protein_synonym              ps
    JOIN protein                      p    ON p.id = ps.protein_id
    JOIN protein_domain_association   pda  ON pda.protein_id = p.id
    JOIN adenylation_domain           ad   ON ad.id = pda.domain_id
    LEFT JOIN substrate_domain_association sda ON sda.domain_id = ad.id
    LEFT JOIN substrate               s    ON s.name = sda.substrate_name
    WHERE ps.synonym = '${t}' COLLATE NOCASE
    ORDER BY pda.domain_number, s.name
    LIMIT 500
              `.replace(/\s+/g, ' ').trim()
              : ''
          );
        }}
        fullWidth
        placeholder="e.g., P48633.1"
      />
    ),
  },
  {
    key: 'species',
    label: 'Substrate specificities by species',
    kind: 'preset',
    defaultQuery: "SELECT DISTINCT species FROM taxonomy LIMIT 500",
    // Editor compiles SQL directly into the main query state
    Editor: ({ presetInput, setPresetInput, setQuery }) => (
      <TextField
        label="Species name"
        value={presetInput}
        onChange={(e) => {
          const v = e.target.value;
          setPresetInput(v);
          const t = v.trim();
          setQuery(
            t
              ? `
    SELECT
      ad.id                    AS domain_id,
      p.id                     AS protein_id,
      t.species                AS species,
      pda.domain_number,
      ad.signature,
      ad.extended_signature,
      s.name                   AS substrate_name,
      s.smiles                 AS substrate_smiles
    FROM taxonomy                     t
    JOIN protein                      p    ON p.taxonomy_id = t.id
    JOIN protein_domain_association   pda  ON pda.protein_id = p.id
    JOIN adenylation_domain           ad   ON ad.id = pda.domain_id
    LEFT JOIN substrate_domain_association sda ON sda.domain_id = ad.id
    LEFT JOIN substrate               s    ON s.name = sda.substrate_name
    WHERE t.species = '${t}' COLLATE NOCASE
    ORDER BY pda.domain_number, s.name
    LIMIT 500
              `.replace(/\s+/g, ' ').trim()
              : ''
          );
        }}
        fullWidth
        placeholder="e.g., Streptomyces coelicolor"
      />
    ),
  },
  {
    key: 'signature',
    label: 'Substrate specificities by A-domain signature (Hamming ≤ N)',
    kind: 'preset',
    defaultQuery: "SELECT * FROM adenylation_domain LIMIT 500",
    // Local state inside editor; compiles SQL into `query` (no parsing back)
    Editor: ({ setQuery }) => {
      const [sig, setSig] = useState('');
      const [maxDist, setMaxDist] = useState(3); // sensible default

      const updateSQL = (s, d) => {
        const clamped = Math.max(0, Math.min(10, Number.isFinite(+d) ? +d : 0));
        const pad = s.toUpperCase().slice(0, 10).padEnd(10, '-'); // GAP token '-'
        // build 10 char-by-char comparisons
        const comps = Array.from({ length: 10 }, (_, i) => {
          const idx = i + 1;
          return `(substr('${pad}',${idx},1) <> substr(UPPER(substr(ad.signature || '----------',1,10)),${idx},1))`;
        }).join(' + ');
        const sql = `
  WITH d AS (
    SELECT
      ad.id,
      ad.signature,
      ad.extended_signature,
      (${comps}) AS hamming
    FROM adenylation_domain ad
  )
  SELECT
    d.id                       AS domain_id,
    p.id                       AS protein_id,
    group_concat(ps.synonym, ', ') AS protein_synonyms
    pda.domain_number,
    ad.signature,
    ad.extended_signature,
    s.name                     AS substrate_name,
    s.smiles                   AS substrate_smiles,
    d.hamming,
  FROM d
  JOIN adenylation_domain           ad  ON ad.id = d.id
  JOIN protein_domain_association   pda ON pda.domain_id = ad.id
  JOIN protein                      p   ON p.id = pda.protein_id
  LEFT JOIN substrate_domain_association sda ON sda.domain_id = ad.id
  LEFT JOIN substrate               s   ON s.name = sda.substrate_name
  LEFT JOIN protein_synonym         ps  ON ps.protein_id = p.id
  WHERE d.hamming <= ${clamped}
  GROUP BY
    d.id, p.id, pda.domain_number,
    ad.signature, ad.extended_signature,
    s.name, s.smiles, d.hamming
  ORDER BY pda.domain_number, s.name
  LIMIT 500
        `.replace(/\s+/g, ' ').trim();
        setQuery(sql);
      };

      return (
        <Stack direction={{ xs: 'column', sm: 'row' }} spacing={1}>
          <TextField
            label="Signature (max 10)"
            value={sig}
            inputProps={{ maxLength: 10 }}
            onChange={(e) => {
              const v = e.target.value;
              setSig(v);
              updateSQL(v, maxDist);
            }}
            fullWidth
            placeholder="e.g., S/T-A-V-I-G-H-D-L"
          />
          <TextField
            label="Max Hamming distance"
            type="number"
            value={maxDist}
            inputProps={{ min: 0, max: 10, step: 1 }}
            onChange={(e) => {
              const v = e.target.value;
              setMaxDist(v);
              updateSQL(sig, v);
            }}
            sx={{ width: 200 }}
          />
        </Stack>
      );
    },
  },
];

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

  const [presetInput, setPresetInput] = useState('');
  const [query, setQuery] = useState(selectedOption.defaultQuery || '');

  // grid state
  const [rows, setRows] = useState([]);
  const [columns, setColumns] = useState([]);
  const [loading, setLoading] = useState(false);
  const [rowCount, setRowCount] = useState(0);
  const [page, setPage] = useState(0);
  const [pageSize, setPageSize] = useState(DEFAULT_PAGE_SIZE);
  const [sortModel, setSortModel] = useState([]);
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

  // Fetch results from the server
  const fetchResults = useCallback(async ({ q, p, ps, sortBy, sortDir }) => {
    if (!q.trim()) {
      toast.warn('Please enter a SQL query.');
      return;
    }
    setLoading(true);
    const reqId = Date.now();
    lastRequestRef.current = reqId;
    try {
      const res = await fetch('/api/sql', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ query: q, page: p, pageSize: ps, sortBy, sortDir }),
      })
      if (!res.ok) throw new Error(await res.text());
      const data = await res.json();
      if (lastRequestRef.current !== reqId) return;
      const cols = data.columns?.length ? data.columns.map((c) => ({ flex: 1, minWidth: 120, ...c })) : buildColumnsFromRows(data.rows || []);
      const withIds = ensureRowIds(data.rows || []);
      setColumns(cols);
      setRows(withIds);
      setRowCount(Number.isFinite(data.total) ? data.total : withIds.length);
    } catch (err) {
      toast.error(`Query failed: ${err.message}`);
      setColumns([]);
      setRows([]);
      setRowCount(0);
    } finally {
      setLoading(false);
    }
  }, [page, selectedOption]);

  // Fetch results when page, pageSize, or activeSort changes
  const onSearch = useCallback(() => {
    setPage(0);
    fetchResults({ q: query, p: 0, ps: pageSize, ...activeSort });
  }, [query, pageSize, activeSort, fetchResults]);

  // Clear all results and reset state
  const clearAll = () => {
    setRows([]);
    setColumns([]);
    setRowCount(0);
    setPage(0);
    setSortModel([]);
    // keep current query text; if preset is active and its editor cleared, query may be ''
  };

  // Export helpers
  const makeDelimited = (cols, data, delim) => {
    const header = cols.map((c) => c.field).join(delim);
    const lines = data.map((r) => cols.map((c) => String(r[c.field] ?? '').replace(/\n|\r/g, ' ')).join(delim));
    return [header, ...lines].join('\n');
  };

  const downloadFile = (content, filename, mime) => {
    const blob = new Blob([content], { type: mime });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    a.remove();
    URL.revokeObjectURL(url);
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
      if (!query.trim()) return toast.warn('Please enter a SQL query first!');
      const res = await fetch('/api/sql', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          query,
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
    if (!query.trim()) return;
    if (rows.length === 0 && rowCount === 0) return;
    fetchResults({ q: query, p: page, ps: pageSize, ...activeSort });
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
              const next = QUERYOPTIONS.find((o) => o.key === nextKey);
              setQuery(next?.defaultQuery || '');
              // clear results on type switch
              setPresetInput('');
              setRows([]);
              setColumns([]);
              setRowCount(0);
              setPage(0);
              setSortModel([]);
            }}
          >
            {QUERYOPTIONS.map((o) => (
              <MenuItem key={o.key} value={o.key}>{o.label}</MenuItem>
            ))}
          </Select>
        </Box>

        {/* Per-option editor:
            - 'free' shows SQL textarea
            - presets render inputs that directly update the main query via setQuery
        */}
        {selectedOption.Editor ? (
          selectedOption.key === 'free'
            ? <selectedOption.Editor query={query} setQuery={setQuery} />
            : <selectedOption.Editor presetInput={presetInput} setPresetInput={setPresetInput} setQuery={setQuery} />
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
            noRowsOverlay: () => (
              <Stack height="100%" alignItems="center" justifyContent="center">
                <Typography variant="body2" color="text.secondary">
                  {rows.length === 0 && !loading ? 'Sorry, no results for query' : ''}
                </Typography>
              </Stack>
            ),
          }}
        />
      </Box>

      {loading && (
        <Box sx={{ position: 'fixed', left: 0, right: 0, top: 0, zIndex: 1200 }}>
          <LinearProgress />
        </Box>
      )}
    </Box>
  );
};

export default QueryDatabase;
