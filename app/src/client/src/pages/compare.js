import React, { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import { useSearchParams } from 'react-router-dom';
import { toast } from 'react-toastify';
import {
    Alert, Autocomplete, Box, Button, Chip, CircularProgress, Divider, IconButton, Menu, MenuItem, Paper,
    TextField, ToggleButton, ToggleButtonGroup, Tooltip, Typography, useTheme,
} from '@mui/material';
import ArrowDownwardIcon from '@mui/icons-material/ArrowDownward';
import ArrowUpwardIcon from '@mui/icons-material/ArrowUpward';
import CloseIcon from '@mui/icons-material/Close';
import DownloadIcon from '@mui/icons-material/Download';
import DragIndicatorIcon from '@mui/icons-material/DragIndicator';
import VerticalAlignTopIcon from '@mui/icons-material/VerticalAlignTop';

import LazyMultiSelect from '../components/LazyMultiSelect';
import { SIMILARITY_RAMP } from '../theme';
import {
    DEFAULT_LABEL_COLUMNS, LABEL_COLUMNS, SIGNATURE_MODES, hammingDistance, moveItem, parseCompareParams,
    residueSimilarity, similarityColors, sortByDistance,
} from '../utils/compare';
import { comparisonToSvg } from '../utils/compareFigure';
import { exportStamp } from '../utils/neighborExport';
import { svgToPng } from '../utils/networkExport';
import { parseSignatureInput } from '../utils/signatures';
import { downloadBlob } from '../utils/zip';

/** Same cap as MAX_ROWS in routes/compare.py. */
const MAX_ROWS = 300;
const CELL_WIDTH = 20;
const CELL_HEIGHT = 24;

/**
 * Fetch reference domains from the server.
 *
 * @param {object} body - {domain_ids, proteins, substrates, species}.
 * @returns {Promise<{rows: Array<object>, truncated: boolean}>} - matching domains.
 */
async function fetchDomains(body) {
    const res = await fetch('/api/compare/domains', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(body),
    });
    const data = await res.json();
    if (!res.ok) throw new Error(data.error || 'failed to load domains');
    return { rows: data.rows.map((r) => ({ ...r, key: `db-${r.id}` })), truncated: data.truncated };
}

let customCounter = 0;
/** A pasted signature as a row; it has an extended signature and nothing else. */
const customRow = (name, extendedSignature) => {
    customCounter += 1;
    return {
        key: `custom-${customCounter}`,
        custom: true,
        name,
        signature: '',
        extended_signature: extendedSignature,
    };
};

/**
 * Compare (extended) signatures of reference domains and pasted signatures,
 * position by position. Each residue is shaded by how physicochemically alike
 * it is to the residue at the same position in the top (reference) row.
 *
 * @returns {React.ReactElement} - The compare page.
 */
const Compare = () => {
    const theme = useTheme();
    const ramp = SIMILARITY_RAMP[theme.palette.mode];
    const [searchParams] = useSearchParams();

    const [rows, setRows] = useState([]);
    const [mode, setMode] = useState('extended');
    const [labelColumns, setLabelColumns] = useState(DEFAULT_LABEL_COLUMNS);
    const [similarity, setSimilarity] = useState(null);
    const [loading, setLoading] = useState(false);

    // pickers for adding database entries
    const [proteins, setProteins] = useState([]);
    const [substrates, setSubstrates] = useState([]);
    const [species, setSpecies] = useState([]);
    const [pasted, setPasted] = useState('');

    // drag and drop: the row being dragged, and where it would land
    const dragFrom = useRef(null);
    const [dropAt, setDropAt] = useState(null);

    useEffect(() => {
        fetch('/api/compare/aa_similarity')
            .then((r) => r.json())
            .then(setSimilarity)
            .catch(() => toast.error('Could not load the amino acid similarity table.'));
    }, []);

    /** Append rows, skipping database domains that are already in the list. */
    // read through a ref so the toasts below stay out of the state updater
    // (StrictMode runs updaters twice, which would show every toast twice)
    const rowsRef = useRef(rows);
    rowsRef.current = rows;
    const addRows = useCallback((incoming, truncated) => {
        const current = rowsRef.current;
        const present = new Set(current.map((r) => r.key));
        const fresh = incoming.filter((r) => !present.has(r.key));
        const room = MAX_ROWS - current.length;
        if (fresh.length > room || truncated) {
            toast.warning(`Only the first ${MAX_ROWS} entries are shown.`);
        }
        if (!fresh.length && incoming.length) toast.info('Those entries are already in the list.');
        const next = [...current, ...fresh.slice(0, Math.max(room, 0))];
        rowsRef.current = next;
        setRows(next);
    }, []);

    // prefill from the URL (e.g. a neighbour list opened from the network page), once
    const prefilled = useRef(false);
    useEffect(() => {
        if (prefilled.current) return;
        prefilled.current = true;
        const { ids, custom } = parseCompareParams(searchParams);
        const customRows = custom.map(({ name, signature }) => customRow(name, signature.toUpperCase()));
        if (!ids.length) {
            if (customRows.length) addRows(customRows, false);
            return;
        }
        setLoading(true);
        fetchDomains({ domain_ids: ids })
            .then(({ rows: dbRows, truncated }) => addRows([...customRows, ...dbRows], truncated))
            .catch((e) => toast.error(e.message))
            .finally(() => setLoading(false));
    }, [searchParams, addRows]);

    const addFromDatabase = async () => {
        setLoading(true);
        try {
            const { rows: dbRows, truncated } = await fetchDomains({ proteins, substrates, species });
            if (!dbRows.length) toast.info('No domains match that selection.');
            addRows(dbRows, truncated);
            setProteins([]);
            setSubstrates([]);
            setSpecies([]);
        } catch (e) {
            toast.error(e.message);
        } finally {
            setLoading(false);
        }
    };

    const addPasted = () => {
        const { queries, errors } = parseSignatureInput(pasted);
        errors.forEach((error) => toast.error(error));
        if (!queries.length) return;
        addRows(queries.map((q, i) => customRow(q.name || `Custom signature ${i + 1}`, q.signature)), false);
        if (!errors.length) setPasted('');
    };

    const field = SIGNATURE_MODES[mode].field;
    const reference = rows[0]?.[field] || '';
    const length = Math.max(0, ...rows.map((r) => (r[field] || '').length));

    // per row: colours per position, and how many positions match the reference
    const rendered = useMemo(() => rows.map((row, rowIndex) => {
        const sequence = row[field] || '';
        let identical = 0;
        const cells = [...sequence].map((residue, position) => {
            const refResidue = reference[position];
            if (residue === refResidue) identical += 1;
            if (rowIndex === 0) return { residue, reference: true };
            const score = residueSimilarity(similarity, residue, refResidue);
            return {
                residue,
                score,
                title: score === null
                    ? `position ${position + 1}: ${residue} (no reference residue to compare)`
                    : `position ${position + 1}: ${residue} vs reference ${refResidue}, similarity ${score.toFixed(2)}`,
                ...(score === null ? {} : similarityColors(score, ramp)),
            };
        });
        return {
            row,
            cells,
            identical,
            compared: Math.min(sequence.length, reference.length),
            distance: hammingDistance(reference, sequence),
        };
    }), [rows, field, reference, similarity, ramp]);

    // download menu anchor
    const [downloadAnchor, setDownloadAnchor] = useState(null);

    /** The table as on screen (order, label columns, mode, shading) as an SVG document. */
    const renderFigure = () => comparisonToSvg({
        columns: labelColumns.map((key) => ({
            label: LABEL_COLUMNS[key].label,
            values: rows.map((row) => LABEL_COLUMNS[key].get(row) || ''),
        })),
        rows: rendered,
        length,
        caption: {
            title: 'Signature comparison',
            subtitle: [
                SIGNATURE_MODES[mode].label,
                `${rows.length} ${rows.length === 1 ? 'entry' : 'entries'}`,
                `reference: ${rows[0]?.name ?? ''}`,
                `shaded by similarity to the reference residue over ${similarity?.n_properties ?? 15} physicochemical descriptors`,
            ].join('; '),
        },
        colors: {
            background: theme.palette.background.paper,
            text: theme.palette.text.primary,
            textSecondary: theme.palette.text.secondary,
            border: theme.palette.divider,
            strongBorder: theme.palette.text.primary,
            referenceFill: theme.palette.surface.sunken,
            rampLow: ramp.low,
            rampHigh: ramp.high,
        },
    });

    const downloadFigure = async (format) => {
        setDownloadAnchor(null);
        const name = `parasect-signature-comparison-${exportStamp()}`;
        try {
            const svg = renderFigure();
            if (format === 'svg') {
                downloadBlob(new Blob([svg], { type: 'image/svg+xml' }), `${name}.svg`);
            } else {
                downloadBlob(await svgToPng(svg), `${name}.png`);
            }
        } catch (e) {
            console.error(e);
            toast.error('Could not create the image.');
        }
    };

    const move = (from, to) => setRows((current) => moveItem(current, from, to));
    const remove = (index) => setRows((current) => current.filter((_, i) => i !== index));

    const onDragOver = (event, index) => {
        event.preventDefault();
        const rect = event.currentTarget.getBoundingClientRect();
        setDropAt(event.clientY < rect.top + rect.height / 2 ? index : index + 1);
    };

    const onDrop = (event) => {
        event.preventDefault();
        const from = dragFrom.current;
        if (from !== null && dropAt !== null) {
            // dropping below itself shifts the target up by the row taken out
            move(from, dropAt > from ? dropAt - 1 : dropAt);
        }
        dragFrom.current = null;
        setDropAt(null);
    };

    const onHandleKey = (event, index) => {
        if (event.key === 'ArrowUp' && index > 0) {
            event.preventDefault();
            move(index, index - 1);
        } else if (event.key === 'ArrowDown' && index < rows.length - 1) {
            event.preventDefault();
            move(index, index + 1);
        }
    };

    const canAdd = proteins.length + substrates.length + species.length > 0;
    const missingInMode = rows.filter((r) => !r[field]).length;

    const headerCell = { px: 1, py: 0.5, textAlign: 'left', fontWeight: 600, fontSize: '0.75rem', whiteSpace: 'nowrap' };
    const labelCell = { px: 1, fontSize: '0.8125rem', whiteSpace: 'nowrap', maxWidth: 240, overflow: 'hidden', textOverflow: 'ellipsis' };

    return (
        <Box sx={{ m: { xs: 2, md: 4 } }}>
            <Typography variant='h4' gutterBottom>Compare signatures</Typography>
            <Typography variant='body1' color='textSecondary' sx={{ maxWidth: 900 }} gutterBottom>
                Put signatures side by side. The top row is the reference: every other residue is shaded by how
                alike it is to the reference residue at the same position, judged on the
                {similarity ? ` ${similarity.n_properties} ` : ' '}
                physicochemical descriptors PARAS featurises signatures with. Signatures all have the same length,
                so no alignment is needed. Drag rows (or focus a handle and use the arrow keys) to reorder them.
            </Typography>

            <Paper variant='outlined' sx={{ p: 2, mt: 2 }}>
                <Typography variant='subtitle1' sx={{ fontWeight: 600, mb: 1 }}>Add from the database</Typography>
                <Box sx={{ display: 'grid', gap: 1.5, gridTemplateColumns: { xs: '1fr', md: '1fr 1fr 1fr auto' }, alignItems: 'start' }}>
                    <LazyMultiSelect field='protein' label='Protein IDs' placeholder='e.g. P48633.1' value={proteins} onChange={setProteins} />
                    <LazyMultiSelect field='substrate' label='Substrates' placeholder='e.g. valine' value={substrates} onChange={setSubstrates} />
                    <LazyMultiSelect field='species' label='Species' placeholder='e.g. Bacillus subtilis' value={species} onChange={setSpecies} />
                    <Button variant='contained' onClick={addFromDatabase} disabled={!canAdd || loading} sx={{ height: 56 }}>
                        Add
                    </Button>
                </Box>
                <Typography variant='subtitle1' sx={{ fontWeight: 600, mt: 2, mb: 1 }}>Or paste extended signatures</Typography>
                <Box sx={{ display: 'flex', gap: 1.5, alignItems: 'stretch', flexWrap: 'wrap' }}>
                    <TextField
                        multiline
                        minRows={2}
                        size='small'
                        value={pasted}
                        onChange={(e) => setPasted(e.target.value)}
                        placeholder={'FASTA, or one per line: name signature'}
                        sx={{ flex: 1, minWidth: 280, '& textarea': { fontFamily: 'monospace' } }}
                    />
                    <Button variant='outlined' onClick={addPasted} disabled={!pasted.trim()}>Add signatures</Button>
                </Box>
            </Paper>

            <Box sx={{ display: 'flex', gap: 2, alignItems: 'center', flexWrap: 'wrap', mt: 3, mb: 1.5 }}>
                <ToggleButtonGroup
                    size='small'
                    exclusive
                    value={mode}
                    onChange={(e, next) => next && setMode(next)}
                    aria-label='Signature to compare'
                >
                    {Object.entries(SIGNATURE_MODES).map(([key, m]) => (
                        <ToggleButton key={key} value={key}>{m.label}</ToggleButton>
                    ))}
                </ToggleButtonGroup>
                <Autocomplete
                    multiple
                    size='small'
                    options={Object.keys(LABEL_COLUMNS)}
                    getOptionLabel={(key) => LABEL_COLUMNS[key].label}
                    value={labelColumns}
                    onChange={(e, next) => setLabelColumns(next)}
                    disableCloseOnSelect
                    renderTags={(selected, getTagProps) => selected.map((key, index) => (
                        <Chip size='small' label={LABEL_COLUMNS[key].label} {...getTagProps({ index })} key={key} />
                    ))}
                    renderInput={(params) => <TextField {...params} label='Label columns (in order)' />}
                    sx={{ minWidth: 320, flex: 1, maxWidth: 640 }}
                />
                <SimilarityLegend ramp={ramp} />
                {rows.length > 2 && (
                    <Box sx={{ display: 'flex', alignItems: 'center', gap: 0.5 }}>
                        <Typography variant='caption' color='textSecondary'>Sort by distance</Typography>
                        <Tooltip title='Closest to the reference first'>
                            <span>
                                <IconButton size='small' aria-label='Sort by distance, closest first' disabled={!reference}
                                    onClick={() => setRows((current) => sortByDistance(current, field, 'asc'))}>
                                    <ArrowUpwardIcon fontSize='small' />
                                </IconButton>
                            </span>
                        </Tooltip>
                        <Tooltip title='Furthest from the reference first'>
                            <span>
                                <IconButton size='small' aria-label='Sort by distance, furthest first' disabled={!reference}
                                    onClick={() => setRows((current) => sortByDistance(current, field, 'desc'))}>
                                    <ArrowDownwardIcon fontSize='small' />
                                </IconButton>
                            </span>
                        </Tooltip>
                    </Box>
                )}
                {rows.length > 0 && (
                    <>
                        <Button
                            size='small'
                            variant='outlined'
                            startIcon={<DownloadIcon fontSize='small' />}
                            onClick={(e) => setDownloadAnchor(e.currentTarget)}
                            aria-haspopup='menu'
                        >
                            Download
                        </Button>
                        <Menu anchorEl={downloadAnchor} open={Boolean(downloadAnchor)} onClose={() => setDownloadAnchor(null)}>
                            <MenuItem onClick={() => downloadFigure('svg')}>SVG (vector)</MenuItem>
                            <MenuItem onClick={() => downloadFigure('png')}>PNG (image)</MenuItem>
                        </Menu>
                        <Button size='small' color='inherit' onClick={() => setRows([])}>Clear all</Button>
                    </>
                )}
            </Box>

            {missingInMode > 0 && (
                <Alert severity='info' sx={{ mb: 1.5 }}>
                    {missingInMode} row{missingInMode === 1 ? ' has' : 's have'} no {SIGNATURE_MODES[mode].label.toLowerCase()}.
                    Pasted signatures are extended only: the 10-residue signature includes a lysine outside the
                    34 extended positions, so it can't be derived from them.
                </Alert>
            )}

            {rows.length > 1 && !reference && (
                <Alert severity='warning' sx={{ mb: 1.5 }}>
                    The reference row ({rows[0].name}) has no {SIGNATURE_MODES[mode].label.toLowerCase()}, so there is
                    nothing to shade against. Move a row that has one to the top.
                </Alert>
            )}

            <Divider />

            {loading && <Box sx={{ py: 3, display: 'flex', justifyContent: 'center' }}><CircularProgress size={28} /></Box>}

            {!loading && rows.length === 0 && (
                <Typography color='textSecondary' sx={{ py: 4 }}>
                    Nothing to compare yet. Add domains from the database or paste signatures above, or open a
                    neighbour list from the network page.
                </Typography>
            )}

            {rows.length > 0 && (
                <Box sx={{ overflowX: 'auto', mt: 1.5 }}>
                    <Box component='table' sx={{ borderCollapse: 'collapse' }} onDragLeave={() => setDropAt(null)}>
                        <thead>
                            <tr>
                                <Box component='th' sx={headerCell} />
                                {labelColumns.map((key) => (
                                    <Box component='th' key={key} sx={headerCell}>{LABEL_COLUMNS[key].label}</Box>
                                ))}
                                <Box component='th' sx={{ ...headerCell, px: 0 }}>
                                    <Box sx={{ display: 'flex' }}>
                                        {Array.from({ length }, (_, i) => (
                                            <Box key={i} sx={{ width: CELL_WIDTH, textAlign: 'center', fontSize: '0.625rem', color: 'text.secondary' }}>
                                                {i === 0 || (i + 1) % 5 === 0 ? i + 1 : ''}
                                            </Box>
                                        ))}
                                    </Box>
                                </Box>
                                <Box component='th' sx={headerCell} title='Hamming distance to the reference: positions that differ'>Distance</Box>
                                <Box component='th' sx={headerCell}>Identical</Box>
                                <Box component='th' sx={headerCell} />
                            </tr>
                        </thead>
                        <tbody>
                            {rendered.map(({ row, cells, identical, compared, distance }, index) => (
                                <Box
                                    component='tr'
                                    key={row.key}
                                    draggable
                                    onDragStart={(e) => {
                                        dragFrom.current = index;
                                        e.dataTransfer.effectAllowed = 'move';
                                        e.dataTransfer.setData('text/plain', row.name);
                                    }}
                                    onDragOver={(e) => onDragOver(e, index)}
                                    onDrop={onDrop}
                                    onDragEnd={() => { dragFrom.current = null; setDropAt(null); }}
                                    sx={{
                                        // the drop position as a line between rows
                                        boxShadow: (t) => (dropAt === index
                                            ? `inset 0 2px 0 ${t.palette.primary.main}`
                                            : dropAt === index + 1 && index === rows.length - 1
                                                ? `inset 0 -2px 0 ${t.palette.primary.main}`
                                                : 'none'),
                                        borderBottom: index === 0 ? '2px solid' : '1px solid',
                                        borderColor: index === 0 ? 'text.primary' : 'divider',
                                        '&:hover .row-actions': { opacity: 1 },
                                    }}
                                >
                                    <Box component='td' sx={{ pl: 0.5, whiteSpace: 'nowrap' }}>
                                        <Tooltip title='Drag to reorder, or focus and use the arrow keys'>
                                            <IconButton
                                                size='small'
                                                aria-label={`Reorder ${row.name}`}
                                                onKeyDown={(e) => onHandleKey(e, index)}
                                                sx={{ cursor: 'grab' }}
                                            >
                                                <DragIndicatorIcon fontSize='small' />
                                            </IconButton>
                                        </Tooltip>
                                        {index === 0 && <Chip size='small' label='Reference' sx={{ ml: 0.5, height: 20, fontSize: '0.6875rem' }} />}
                                    </Box>
                                    {labelColumns.map((key) => {
                                        const value = LABEL_COLUMNS[key].get(row) || (row.custom && key === 'name' ? row.name : '');
                                        return (
                                            <Box component='td' key={key} sx={labelCell} title={value}>
                                                {value || <Box component='span' sx={{ color: 'text.disabled' }}>–</Box>}
                                            </Box>
                                        );
                                    })}
                                    <Box component='td' sx={{ px: 0, py: 0.5 }}>
                                        {cells.length ? (
                                            <Box sx={{ display: 'flex', fontFamily: 'monospace', fontSize: '0.8125rem' }}>
                                                {cells.map((cell, position) => (
                                                    <Box
                                                        // eslint-disable-next-line react/no-array-index-key
                                                        key={position}
                                                        title={cell.title}
                                                        sx={{
                                                            width: CELL_WIDTH,
                                                            height: CELL_HEIGHT,
                                                            lineHeight: `${CELL_HEIGHT}px`,
                                                            textAlign: 'center',
                                                            backgroundColor: cell.reference ? 'surface.sunken' : cell.fill,
                                                            color: cell.reference ? 'text.primary' : cell.text,
                                                            fontWeight: cell.reference ? 700 : 500,
                                                        }}
                                                    >
                                                        {cell.residue}
                                                    </Box>
                                                ))}
                                            </Box>
                                        ) : (
                                            <Typography variant='caption' color='textSecondary' sx={{ px: 1 }}>not available</Typography>
                                        )}
                                    </Box>
                                    <Box component='td' sx={{ ...labelCell, fontVariantNumeric: 'tabular-nums', textAlign: 'right' }}>
                                        {index === 0 || distance === null ? '' : distance}
                                    </Box>
                                    <Box component='td' sx={{ ...labelCell, fontVariantNumeric: 'tabular-nums', color: 'text.secondary' }}>
                                        {index === 0 || !compared ? '' : `${identical}/${compared}`}
                                    </Box>
                                    <Box component='td' sx={{ whiteSpace: 'nowrap', pr: 0.5 }}>
                                        <Box className='row-actions' sx={{ opacity: { xs: 1, md: 0.35 }, transition: 'opacity 120ms', '&:focus-within': { opacity: 1 } }}>
                                            {index > 0 && (
                                                <Tooltip title='Make this the reference'>
                                                    <IconButton size='small' aria-label={`Make ${row.name} the reference`} onClick={() => move(index, 0)}>
                                                        <VerticalAlignTopIcon fontSize='small' />
                                                    </IconButton>
                                                </Tooltip>
                                            )}
                                            <Tooltip title='Remove'>
                                                <IconButton size='small' aria-label={`Remove ${row.name}`} onClick={() => remove(index)}>
                                                    <CloseIcon fontSize='small' />
                                                </IconButton>
                                            </Tooltip>
                                        </Box>
                                    </Box>
                                </Box>
                            ))}
                        </tbody>
                    </Box>
                </Box>
            )}
        </Box>
    );
};

/** Gradient key for the residue shading. */
const SimilarityLegend = ({ ramp }) => (
    <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
        <Typography variant='caption' color='textSecondary'>Dissimilar</Typography>
        <Box
            sx={{
                width: 120,
                height: 12,
                borderRadius: 0.5,
                border: '1px solid',
                borderColor: 'divider',
                background: `linear-gradient(to right, ${ramp.low}, ${ramp.high})`,
            }}
        />
        <Typography variant='caption' color='textSecondary'>Identical</Typography>
    </Box>
);

export default Compare;
