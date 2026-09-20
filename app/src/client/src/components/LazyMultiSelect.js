import React, { useEffect, useMemo, useRef, useState } from 'react';
import { Autocomplete, Box, Chip, CircularProgress, TextField, Typography } from '@mui/material';

/** How many options to ask the server for at a time. */
const FETCH_LIMIT = 100;

/**
 * Multi-select filter backed by /api/sql/options.
 *
 * The values are looked up as you type rather than pulled down in full: there
 * are already ~1,900 distinct protein identifiers in the reference database,
 * and mounting that many rows in a listbox is what made the old dropdowns
 * heavy. The server returns the most-used matches first and caps the list.
 *
 * @param {object} props - component props.
 * @param {'substrate'|'protein'|'species'} props.field - which option list to search.
 * @param {string[]} props.value - currently selected values.
 * @param {(next: string[]) => void} props.onChange - called with the new selection.
 * @param {string} props.label - field label.
 * @param {string} props.placeholder - input placeholder.
 * @returns {React.ReactElement} - the filter control.
 */
const LazyMultiSelect = ({ field, value, onChange, label, placeholder }) => {
    const [input, setInput] = useState('');
    const [options, setOptions] = useState([]);
    const [loading, setLoading] = useState(false);
    const [truncated, setTruncated] = useState(false);
    const requestRef = useRef(0);

    useEffect(() => {
        const query = input.trim();
        setLoading(true);
        const handle = setTimeout(() => {
            const reqId = requestRef.current + 1;
            requestRef.current = reqId;
            fetch(`/api/sql/options?field=${field}&q=${encodeURIComponent(query)}&limit=${FETCH_LIMIT}`)
                .then((r) => {
                    if (!r.ok) throw new Error('option lookup failed');
                    return r.json();
                })
                .then((data) => {
                    // ignore anything but the newest request, so a slow early
                    // response can't overwrite the results for what's typed now
                    if (requestRef.current !== reqId) return;
                    const received = data.options || [];
                    setOptions(received);
                    setTruncated(received.length >= FETCH_LIMIT);
                })
                .catch(() => {
                    if (requestRef.current === reqId) setOptions([]);
                })
                .finally(() => {
                    if (requestRef.current === reqId) setLoading(false);
                });
        }, 250);
        return () => clearTimeout(handle);
    }, [field, input]);

    const countByValue = useMemo(
        () => new Map(options.map((o) => [o.value, o.count])),
        [options],
    );

    // keep already-selected values in the list so their chips survive a re-search
    const selectable = useMemo(() => {
        const names = options.map((o) => o.value);
        const missing = value.filter((v) => !names.includes(v));
        return [...missing, ...names];
    }, [options, value]);

    return (
        <Box>
            <Autocomplete
                multiple
                disableCloseOnSelect
                filterSelectedOptions
                options={selectable}
                value={value}
                inputValue={input}
                onInputChange={(e, next, reason) => { if (reason !== 'reset') setInput(next); }}
                onChange={(e, next) => onChange(next)}
                // the server already searched; re-filtering here would hide
                // matches it found by a different rule
                filterOptions={(opts) => opts}
                loading={loading}
                noOptionsText={loading ? 'Searching...' : 'No matches'}
                renderOption={(props, option) => (
                    <li {...props} key={option}>
                        <Box component='span' sx={{ flex: 1, overflow: 'hidden', textOverflow: 'ellipsis' }}>
                            {option}
                        </Box>
                        {countByValue.has(option) && (
                            <Box component='span' sx={{ ml: 1, color: 'text.secondary', fontSize: '0.75rem' }}>
                                {countByValue.get(option)}
                            </Box>
                        )}
                    </li>
                )}
                renderTags={(selected, getTagProps) => selected.map((option, index) => (
                    <Chip size='small' label={option} {...getTagProps({ index })} key={option} />
                ))}
                renderInput={(inputParams) => (
                    <TextField
                        {...inputParams}
                        label={label}
                        placeholder={value.length ? '' : placeholder}
                        InputProps={{
                            ...inputParams.InputProps,
                            endAdornment: (
                                <>
                                    {loading ? <CircularProgress size={18} /> : null}
                                    {inputParams.InputProps.endAdornment}
                                </>
                            ),
                        }}
                    />
                )}
                fullWidth
            />
            {truncated && (
                <Typography variant='caption' color='textSecondary' sx={{ mt: 0.5, display: 'block' }}>
                    Showing the {FETCH_LIMIT} most-used matches - keep typing to narrow.
                </Typography>
            )}
        </Box>
    );
};

export default LazyMultiSelect;
