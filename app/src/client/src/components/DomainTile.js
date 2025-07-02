import React, { useState, useEffect } from 'react';
import { toast } from 'react-toastify';
import {
    Box,
    Button,
    MenuItem,
    Select,
    FormControl,
    Autocomplete,
    TextField,
    FormControlLabel,
    Checkbox,
    Popper,
} from '@mui/material';
import { FaFingerprint, FaCopy } from 'react-icons/fa';

import SmileDrawerContainer from './SmilesDrawer';

/**
 * Component to display the results of the prediction.
 *
 * @param {Object} props - The component props.
 * @param {Object} props.result - The result object.
 * @returns {React.ReactElement} - The result tile component.
 */
const DomainTile = ({ result }) => {
    const [selectedPrediction] = useState(result['predictions'][0]);
    const [selectedSubstrate, setSelectedSubstrate] = useState(null);

    const [useCustomSmiles, setUseCustomSmiles] = useState(false);
    const [customSmiles, setCustomSmiles] = useState('');
    const [newSubstrateName, setNewSubstrateName] = useState(null)

    const [smilesOptions, setSmilesOptions] = useState([]);

    {/* Load SMILES for substrate assignment */}
    useEffect(() => {
        fetch('/smiles.tsv')
            .then((res) => res.text())
            .then((text) => {
                const lines = text.trim().split('\n');
                const [headerLine, ...dataLines] = lines;
                const headers = ["substrate_name", "substrate_smiles"]

                const data = dataLines.map(line => {
                    const values = line.split('\t');
                    return headers.reduce((obj, key, i) => {
                        obj[key] = values[i];
                        return obj;
                    }, {});
                });

                setSmilesOptions(data); // [{ substrate, smiles }, ...]
            })

            .catch((err) => {
                console.error('Failed to load smiles.tsv:', err);
            });

    }, []);

    // Sort the smilesOptions by PARAS prediction; stick all SMILES not in the model at the end:
    const [sortedOptions, setSortedOptions] = useState([]);

    useEffect(() => {
          if (!smilesOptions.length || !result["predictions"]) return;

          const preferredOrder = result["predictions"].map(p => p["substrate_name"]);
          const orderMap = new Map(preferredOrder.map((name, i) => [name, i]));

          const sorted = [...smilesOptions].sort((a, b) => {
            const aIndex = orderMap.has(a["substrate_name"]) ? orderMap.get(a["substrate_name"]) : Infinity;
            const bIndex = orderMap.has(b["substrate_name"]) ? orderMap.get(b["substrate_name"]) : Infinity;
            if (aIndex !== bIndex) return aIndex - bIndex;
            return a["substrate_name"].localeCompare(b["substrate_name"]);
          });

          setSortedOptions(sorted);
        }, [smilesOptions, result["predictions"]]);

    return (
        <Box
            sx={{
                minWidth: '500px',
                maxWidth: '500px',
                borderRadius: '11px',
                boxShadow: '0px 4px 10px rgba(100, 84, 31, 0.5)',
            }}
        >
            {/* header */}
            <Box
                sx={{
                    backgroundColor: 'secondary.main',
                    color: 'black.main',
                    padding: '14px 8px',
                    display: 'flex',
                    marginBottom: 1,
                    borderTopLeftRadius: '10px',
                    borderTopRightRadius: '10px',
                    justifyContent: 'center',
                }}
            >
                {`Domain ${result['domain_nr']} (${result['domain_start']}-${result['domain_end']})`}
            </Box>

            {/* domain signature */}
            {result['domain_signature'].length > 0 && (
                <Box
                    sx={{
                        padding: 1,
                        display: 'flex',
                        justifyContent: 'left',
                    }}
                >
                    <FaFingerprint size="1.5em" style={{ marginRight: '5px' }} />
                    {result['domain_signature']}
                </Box>
            )}

            {/* spacing */}
            <Box sx={{ height: 16 }} />

            {/* prediction + substrate selector + visualizer */}
            <Box
                sx={{
                    display: 'flex',
                    alignItems: 'flex-start',
                    paddingX: 1,
                    gap: 2,
                }}
            >
                <Box sx={{ flex: 1 }}>
                    <Box sx={{ display: 'flex', flexDirection: 'column', gap: 2 }}>
                        {/* PARAS prediction */}
                        <Box>
                            <Box sx={{ fontWeight: 500, marginBottom: 0.5 }}>
                                PARAS prediction
                            </Box>
                            <Box
                                sx={{
                                    border: '1px solid #ccc',
                                    borderRadius: '4px',
                                    padding: '8px',
                                    minHeight: '40px',
                                    display: 'flex',
                                    alignItems: 'center',
                                    justifyContent: 'flex-start',
                                    textAlign: 'left',
                                }}
                            >
                                {selectedPrediction['substrate_name']} (
                                {selectedPrediction['probability'].toFixed(2)})
                            </Box>
                        </Box>

                        {/* Substrate */}
                        <Box>
                            <Box sx={{ fontWeight: 500, marginBottom: 0.5 }}>
                                Substrate 1
                            </Box>
                            <Autocomplete
                                options={sortedOptions}
                                getOptionLabel={(option) => option['substrate_name']}
                                renderInput={(params) => (
                                    <TextField {...params} variant="outlined" />
                                )}
                                value={selectedSubstrate || null}
                                onChange={(event, newValue) => {
                                    setSelectedSubstrate(newValue);
                                }}
                                isOptionEqualToValue={(option, value) =>
                                    option['substrate_name'] === value['substrate_name']
                                }
                                clearOnEscape
                                sx={{ borderRadius: '0' }}
                                PopperComponent={(props) => (
                                    <Popper {...props}
                                            placement="bottom-start"
                                            modifiers={[{
                                                name: 'flip',
                                                enabled: false
                                            }]}
                                    />
                                )}
                            />
                        </Box>

                        <Box>
                            <FormControlLabel
                                control={
                                    <Checkbox
                                        checked={useCustomSmiles}
                                        onChange={(e) => {
                                            setUseCustomSmiles(e.target.checked);
                                            setCustomSmiles('');
                                        }}
                                    />
                                }
                                label="New substrate"
                            />
                            {useCustomSmiles && (
                                <Box sx={{ mt: 1 }}>
                                    <TextField
                                        label="Enter substrate name"
                                        fullWidth
                                        value={newSubstrateName}
                                        onChange={(e) => {
                                            const value = e.target.value;
                                            setNewSubstrateName(value);

                                        }}

                                    />
                                </Box>
                            )}

                            {useCustomSmiles && (
                                <Box sx={{ mt: 1 }}>
                                    <TextField
                                        label="Enter SMILES"
                                        fullWidth
                                        value={customSmiles}
                                        onChange={(e) => {
                                            const value = e.target.value;
                                            setCustomSmiles(value);

                                        }}

                                    />
                                </Box>
                            )}

                        </Box>
                    </Box>
                </Box>

                {/* substrate visualizer */}
                {useCustomSmiles && customSmiles ? (
                    <SmileDrawerContainer
                        identifier={`${result['domain_name']}-${result['domain_nr']}`}
                        smilesStr={customSmiles}
                        height={200}
                        width={200}
                    />
                ) : selectedSubstrate ? (
                    <SmileDrawerContainer
                        identifier={`${result['domain_name']}-${result['domain_nr']}`}
                        smilesStr={selectedSubstrate['substrate_smiles']}
                        height={200}
                        width={200}
                    />
                ) : (
                    <Box
                        sx={{
                            height: 200,
                            width: 200,
                            display: 'flex',
                            alignItems: 'center',
                            justifyContent: 'center',
                            color: 'gray',
                            fontStyle: 'italic',
                        }}
                    >
                        No substrate selected
                    </Box>
                )}
            </Box>

            {/* Copy buttons section */}
            <Box>
                <Box
                    sx={{
                        display: 'flex',
                        justifyContent: 'center',
                        alignItems: 'center',
                        gap: '1px',
                    }}
                >
                    <Button
                        variant="contained"
                        color="primary"
                        onClick={() => {
                            navigator.clipboard.writeText(result['domain_sequence']);
                            toast.success('Copied the domain amino acid sequence to clipboard!');
                        }}
                        sx={{
                            flexGrow: 1,
                            width: '50%',
                            borderRadius: '0',
                            borderBottom: '1px solid white',
                        }}
                        disabled={result['domain_sequence'].length === 0}
                    >
                        <FaCopy style={{ marginRight: '5px', fill: 'white' }} />
                        Sequence
                    </Button>
                    <Button
                        variant="contained"
                        color="primary"
                        onClick={() => {
                            navigator.clipboard.writeText(result['domain_signature']);
                            toast.success('Copied the domain signature to clipboard!');
                        }}
                        sx={{
                            flexGrow: 1,
                            width: '50%',
                            borderRadius: '0',
                            borderBottom: '1px solid white',
                        }}
                        disabled={result['domain_signature'].length === 0}
                    >
                        <FaCopy style={{ marginRight: '5px', fill: 'white' }} />
                        Signature
                    </Button>
                </Box>

                <Box
                    sx={{
                        display: 'flex',
                        justifyContent: 'center',
                        alignItems: 'center',
                        gap: '1px',
                    }}
                >
                    <Button
                        variant="contained"
                        color="primary"
                        onClick={() => {
                            navigator.clipboard.writeText(selectedPrediction['substrate_smiles']);
                            toast.success('Copied the substrate SMILES to clipboard!');
                        }}
                        sx={{
                            flexGrow: 1,
                            width: '50%',
                            borderRadius: '0',
                            borderBottomLeftRadius: '10px',
                        }}
                        disabled={selectedPrediction['substrate_smiles'].length === 0}
                    >
                        <FaCopy style={{ marginRight: '5px', fill: 'white' }} />
                        SMILES
                    </Button>
                    <Button
                        variant="contained"
                        color="primary"
                        onClick={() => {
                            navigator.clipboard.writeText(result['domain_extended_signature']);
                            toast.success('Copied the domain extended signature to clipboard!');
                        }}
                        sx={{
                            flexGrow: 1,
                            width: '50%',
                            borderRadius: '0',
                            borderBottomRightRadius: '10px',
                        }}
                        disabled={result['domain_extended_signature'].length === 0}
                    >
                        <FaCopy style={{ marginRight: '5px', fill: 'white' }} />
                        Ext. signature
                    </Button>
                </Box>
            </Box>
        </Box>
    );
};

export default DomainTile;