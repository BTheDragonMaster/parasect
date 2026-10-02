import React, { useState } from 'react';
import { Alert, Box, Button, Chip, CircularProgress, TextField, Tooltip, Typography } from '@mui/material';
import DownloadIcon from '@mui/icons-material/Download';

import { MAX_QUERIES, describePlacement, parseSignatureInput } from '../utils/signatures';

/** Parse problems shown at once; the rest are counted, not listed. */
const ERRORS_SHOWN = 5;

/**
 * A dot with a ring around it: the mark a placed signature gets in the graph.
 *
 * @param {object} props - component props.
 * @param {string} props.color - the query mark colour.
 * @returns {React.ReactElement} - the swatch.
 */
export const QuerySwatch = ({ color }) => (
    <Box
        component='span'
        sx={{
            display: 'inline-block',
            width: 9,
            height: 9,
            borderRadius: '50%',
            flexShrink: 0,
            backgroundColor: color,
            boxShadow: (theme) => `0 0 0 1.5px ${theme.palette.background.paper}, 0 0 0 3px ${color}`,
            ml: '7px !important',
        }}
    />
);

/**
 * Sidebar section of the network page where users place their own extended
 * signatures: paste them, see one chip per placed signature, pick one to list
 * its neighbours, and download everything.
 *
 * @param {object} props - component props.
 * @param {Array<{key: string, name: string, signature: string}>} props.queries - placed signatures.
 * @param {Object<string, object>} props.placements - server placement per query key, at the current threshold.
 * @param {string|null} props.selectedKey - the signature whose neighbours are listed.
 * @param {number} props.threshold - the current threshold.
 * @param {boolean} props.placing - whether a placement request is in flight.
 * @param {string|null} props.error - the last placement error.
 * @param {boolean} props.downloading - whether the results zip is being built.
 * @param {string} props.queryColor - the query mark colour.
 * @param {(parsed: Array<{name: string|null, signature: string}>) => Promise<boolean>} props.onPlace -
 *     place new signatures; resolves true once they are in the graph.
 * @param {(key: string) => void} props.onSelect - list a signature's neighbours and pan to it.
 * @param {(key: string) => void} props.onRemove - take one signature out.
 * @param {() => void} props.onClear - take them all out.
 * @param {() => void} props.onDownload - download tables and figure for every placed signature.
 * @returns {React.ReactElement} - the panel.
 */
const SignatureQueryPanel = ({
    queries, placements, selectedKey, threshold, placing, error, downloading, queryColor,
    onPlace, onSelect, onRemove, onClear, onDownload,
}) => {
    const [text, setText] = useState('');
    const [parseErrors, setParseErrors] = useState([]);

    const submit = async () => {
        const { queries: parsed, errors } = parseSignatureInput(text);
        if (!parsed.length && !errors.length) {
            setParseErrors(['Paste at least one signature first.']);
            return;
        }
        setParseErrors(errors);
        // all or nothing: placing the good half of a paste quietly would leave
        // the bad half looking placed too
        if (errors.length) return;
        if (await onPlace(parsed)) setText('');
    };

    const placedCount = queries.filter((q) => placements[q.key]).length;

    return (
        <Box>
            <Box sx={{ display: 'flex', alignItems: 'baseline', justifyContent: 'space-between' }}>
                <Typography variant='subtitle2'>Your signatures</Typography>
                {queries.length > 0 && (
                    <Typography variant='caption' color='textSecondary'>
                        {queries.length} of {MAX_QUERIES}
                    </Typography>
                )}
            </Box>
            <Typography variant='caption' color='textSecondary' display='block' sx={{ mb: 1 }}>
                Paste 34-residue extended signatures, as FASTA or one per line with an optional name in front,
                to see which clusters they would join. They are placed at the current threshold and not stored
                anywhere.
            </Typography>

            <TextField
                size='small'
                fullWidth
                multiline
                minRows={3}
                maxRows={8}
                placeholder={'>my_domain\nLDASFDASLFEVWGALTAGA...'}
                value={text}
                onChange={(e) => { setText(e.target.value); setParseErrors([]); }}
                onKeyDown={(e) => {
                    if (e.key === 'Enter' && (e.metaKey || e.ctrlKey)) {
                        e.preventDefault();
                        submit();
                    }
                }}
                inputProps={{ spellCheck: false, 'aria-label': 'Extended signatures to place' }}
                InputProps={{ sx: { fontFamily: 'monospace', fontSize: '0.8rem' } }}
                sx={{ mb: 1 }}
            />

            {parseErrors.length > 0 && (
                <Alert severity='warning' variant='outlined' sx={{ py: 0, mb: 1, fontSize: '0.75rem' }}>
                    {parseErrors.slice(0, ERRORS_SHOWN).map((message) => (
                        <Box key={message} sx={{ overflowWrap: 'anywhere' }}>{message}</Box>
                    ))}
                    {parseErrors.length > ERRORS_SHOWN && (
                        <Box>and {parseErrors.length - ERRORS_SHOWN} more.</Box>
                    )}
                </Alert>
            )}
            {error && (
                <Alert severity='error' variant='outlined' sx={{ py: 0, mb: 1, fontSize: '0.75rem' }}>
                    {error}
                </Alert>
            )}

            <Box sx={{ display: 'flex', gap: 1, mb: 1 }}>
                <Button
                    size='small'
                    variant='contained'
                    fullWidth
                    onClick={submit}
                    disabled={placing || !text.trim()}
                    startIcon={placing ? <CircularProgress size={14} color='inherit' /> : null}
                >
                    {placing ? 'Placing...' : 'Place in network'}
                </Button>
                <Button
                    size='small'
                    variant='outlined'
                    fullWidth
                    onClick={onDownload}
                    disabled={!placedCount || downloading || placing}
                    startIcon={downloading
                        ? <CircularProgress size={14} color='inherit' />
                        : <DownloadIcon fontSize='small' />}
                >
                    {downloading ? 'Packing...' : 'Results'}
                </Button>
            </Box>

            {queries.length > 0 && (
                <>
                    <Box sx={{ display: 'flex', flexWrap: 'wrap', gap: 0.5, mb: 0.5 }}>
                        {queries.map((q) => {
                            const isSelected = q.key === selectedKey;
                            return (
                                <Tooltip key={q.key} title={describePlacement(placements[q.key], threshold)}>
                                    <Chip
                                        size='small'
                                        variant='outlined'
                                        icon={<QuerySwatch color={queryColor} />}
                                        label={q.name}
                                        onClick={() => onSelect(q.key)}
                                        onDelete={() => onRemove(q.key)}
                                        sx={{
                                            maxWidth: '100%',
                                            color: 'text.primary',
                                            fontWeight: isSelected ? 700 : 400,
                                            backgroundColor: isSelected ? 'action.selected' : 'transparent',
                                            opacity: placements[q.key] ? 1 : 0.6,
                                        }}
                                    />
                                </Tooltip>
                            );
                        })}
                    </Box>
                    <Button size='small' onClick={onClear} disabled={placing}>
                        Remove all
                    </Button>
                </>
            )}
        </Box>
    );
};

export default SignatureQueryPanel;
