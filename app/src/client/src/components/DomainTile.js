import React, { useState, useEffect, useMemo } from 'react';
import {
    ExpandMore,
    ExpandLess
} from '@mui/icons-material';

import {
    Box,
    Button,
    Autocomplete,
    TextField,
    FormControlLabel,
    Checkbox,
    Popper,
    Collapse,
    IconButton,
    Typography,
} from '@mui/material';
import { FaFingerprint } from 'react-icons/fa';

import SmileDrawerContainer from './SmilesDrawer';
import SmilesChecker from './SmilesChecker';

/**
 * Component to display the results of the prediction.
 *
 * @param {Object} props - The component props.
 * @param {Object} props.result - The result object.
 * @returns {React.ReactElement} - The result tile component.
 */
const DomainTile = ({ result , domainIndex, protein_name, onAnnotationChange}) => {
    const [expanded, setExpanded] = useState(false);  // Start expanded
    const toggleExpanded = () => setExpanded(prev => !prev);
    const [nameWarnings, setNameWarnings] = useState({});
    const [nameValidity, setNameValidity] = useState({});

    const parasResult = result["paras_result"];
    const sequence = parasResult["domain_sequence"];
    const signature = parasResult["domain_signature"];
    const extendedSignature = parasResult["domain_extended_signature"];

    const sequenceMatches = result["sequence_matches"];
    const hasSequenceMatch = sequenceMatches && sequenceMatches.length > 0;
    const firstSequenceMatch = hasSequenceMatch ? sequenceMatches[0] : null;
    const matchSynonym = hasSequenceMatch ? firstSequenceMatch["synonyms"][0]["synonym"] : null;
    const matchedSubstrates = hasSequenceMatch ? firstSequenceMatch["substrates"] : null;
    const [overrideSubstrates, setOverrideSubstrates] = useState(false);


    const [proteinName, setProteinName] = useState(protein_name || '');
    const [isDuplicateDomain, setIsDuplicateDomain] = useState(false);
    const [annotatedSubstrates, setAnnotatedSubstrates] = useState([]);

    const domainName = useMemo(() => {
        return matchSynonym != null
            ? matchSynonym
            : `${protein_name}.A${domainIndex}`;
    }, [matchSynonym, protein_name, domainIndex]);

    const domainSynonym = useMemo( () => {
        return `${proteinName}.A${domainIndex}`;}, [protein_name, domainIndex]);

    useEffect(() => {
    const hasValidSubstrate = substrates?.some(
        sub => sub.substrateName && sub.substrateSmiles
    );

    if (hasValidSubstrate) {
        updateSubstrateField(substrates);
    }
}, [isDuplicateDomain]);

    useEffect(() => {
        const fetchDuplicateStatus = async () => {
            if (domainName) {
                try {
                    const response = await fetch('/api/check_domain_name', {
                        method: 'POST',
                        headers: {
                            'Content-Type': 'application/json',
                        },
                        body: JSON.stringify({ domain_name: domainName }),
                    });

                    const data = await response.json();
                    setIsDuplicateDomain(data.domain_in_dataset ?? false);
                } catch (error) {
                    console.error('Error checking domain name:', error);
                    setIsDuplicateDomain(false);
                }
            }
        };

        fetchDuplicateStatus();
    }, [domainName]);

    {/* Initialize with one substrate entry, default to no selection & no custom smiles */
    }
    const [substrates, setSubstrates] = useState([
        {
            selectedSubstrate: null,
            useCustomSmiles: false,
            customSmiles: '',
            newSubstrateName: '',
            substrateName: '',
            substrateSmiles: '',
            annotationType: null,
            sequence: sequence,
            signature: signature,
            extendedSignature: extendedSignature},
    ]);

    const [smilesOptions, setSmilesOptions] = useState([]);

    // Helper to compare known vs selected (dropdown-based) substrates
    const matchedNames = new Set((matchedSubstrates || []).map(s => s.name));

    const selectedDropdownNames = new Set(
        substrates
            .filter(sub => !sub.useCustomSmiles && sub.selectedSubstrate)
            .map(sub => sub.selectedSubstrate.name)
    );

    const allSelectedAreMatched = selectedDropdownNames.size === matchedNames.size &&
        [...selectedDropdownNames].every(name => matchedNames.has(name));

    const showNoUpdateMessage = hasSequenceMatch && allSelectedAreMatched;

    const recalculateAnnotationTypes = (substrates) => {
    const matchedNames = new Set((matchedSubstrates || []).map(s => s.name));
    const selectedDropdownNames = new Set(
        substrates
            .filter(sub => !sub.useCustomSmiles && sub.selectedSubstrate)
            .map(sub => sub.selectedSubstrate.name)
    );

    const allMatch =
        selectedDropdownNames.size === matchedNames.size &&
        [...selectedDropdownNames].every(name => matchedNames.has(name));

    return substrates.map((sub) => {
        let annotationType = null;

        // 🚨 Overwrite logic for duplicate domains without sequence match
        if (!hasSequenceMatch && isDuplicateDomain) {
            annotationType = 'duplicate_entry';
        } else if (!sub.useCustomSmiles && sub.selectedSubstrate) {
            if (hasSequenceMatch && allMatch) {
                annotationType = 'no_update';
            } else if (hasSequenceMatch) {
                annotationType = 'correction';
            } else {
                annotationType = 'new_entry';
            }
        } else if (sub.useCustomSmiles) {
            annotationType = hasSequenceMatch ? 'correction' : 'new_entry';
        }

        return { ...sub, annotationType };
    });
};

    {/* Load SMILES for substrate assignment */
    }
    useEffect(() => {
          fetch('/api/get_substrates')
            .then(res => res.json())
            .then(data => {
              setSmilesOptions(data);
            })
            .catch(err => {
              console.error('Failed to load smiles data from API:', err);
            });
        }, []);

    {/* Sort smilesOptions by PARAS predictions order */
    }
    const [sortedOptions, setSortedOptions] = useState([]);

    useEffect(() => {
        if (!smilesOptions.length || !parasResult["predictions"]) return;

        const preferredOrder = parasResult["predictions"].map(p => p["substrate_name"]);
        const orderMap = new Map(preferredOrder.map((name, i) => [name, i]));

        const sorted = [...smilesOptions].sort((a, b) => {
            const aIndex = orderMap.has(a["name"]) ? orderMap.get(a["name"]) : Infinity;
            const bIndex = orderMap.has(b["name"]) ? orderMap.get(b["name"]) : Infinity;
            if (aIndex !== bIndex) return aIndex - bIndex;
            return a["name"].localeCompare(b["name"]);
        });

        setSortedOptions(sorted);
    }, [smilesOptions, parasResult["predictions"]]);

    {/* The first prediction for showing */}
    const selectedPrediction = parasResult['predictions'][0];



    {/* Handlers to update substrate entries */}
    const updateSubstrateField = (index, field, value) => {
    setSubstrates((prev) => {
        const newSubs = [...prev];
        const updated = { ...newSubs[index], [field]: value };

        if (field === 'useCustomSmiles') {
            if (value === true) {
                updated.selectedSubstrate = null;
                updated.substrateName = '';
                updated.substrateSmiles = '';
            } else {

                updated.customSmiles = '';
                updated.newSubstrateName = '';
                updated.substrateName = '';
                updated.substrateSmiles = '';
            }
        }

        if (field === 'selectedSubstrate' && value !== null) {
            updated.useCustomSmiles = false;
            updated.customSmiles = '';
            updated.newSubstrateName = '';
            updated.substrateName = value.name;
            updated.substrateSmiles = value.smiles;
        }

        newSubs[index] = updated;

        const finalSubs = recalculateAnnotationTypes(newSubs);
        onAnnotationChange?.(domainName, finalSubs);
        return finalSubs;
    });
};
    useEffect(() => {
    substrates.forEach((sub, i) => {
        const name = sub.newSubstrateName?.trim();
        const smiles = sub.customSmiles?.trim();

        if (sub.useCustomSmiles && name) {
            fetch('/api/check_substrate_name', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ substrate_name: name }),
            })
                .then((res) => res.json())
                .then((data) => {
                    const duplicate = Array.isArray(data) && data.length > 0;

                    setNameWarnings((prev) => ({
                        ...prev,
                        [i]: duplicate
                            ? 'Substrate name already exists in the dataset.'
                            : null,
                    }));

                    setNameValidity((prev) => ({
                        ...prev,
                        [i]: !duplicate,
                    }));

                    setSubstrates((prev) => {
                        const updated = [...prev];
                        const item = { ...updated[i] };

                        if (!duplicate && name && smiles) {
                            item.substrateName = name;
                            item.substrateSmiles = smiles;
                        } else {
                            item.substrateName = '';
                            item.substrateSmiles = '';
                        }

                        updated[i] = item;
                        onAnnotationChange?.(domainName, updated);
                        return updated;
                    });
                })
                .catch((err) => {
                    console.error('Failed to check substrate name:', err);
                    setNameWarnings((prev) => ({
                        ...prev,
                        [i]: 'Error checking substrate name.',
                    }));
                });
        }
    });
}, [substrates.map((s) => `${s.newSubstrateName}|${s.customSmiles}`).join('|')]);

    useEffect(() => {
    substrates.forEach((sub, i) => {
        const name = sub.newSubstrateName?.trim();
        if (sub.useCustomSmiles && name) {
            fetch('/api/check_substrate_name', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ substrate_name: name }),
                })
                .then((res) => res.json())
                .then((data) => {
                    const isDuplicate = Array.isArray(data) && data.length > 0;
                    setNameWarnings((prev) => ({
                        ...prev,
                        [i]: isDuplicate ? `Substrate name already exists in the dataset.` : null,
                    }));
                    setNameValidity((prev) => ({
                        ...prev,
                        [i]: !isDuplicate,
                    }));
                })
                .catch((err) => {
                    console.error('Failed to check substrate name:', err);
                    setNameWarnings((prev) => ({
                        ...prev,
                        [i]: 'Error checking substrate name.',
                    }));
                    setNameValidity((prev) => ({
                        ...prev,
                        [i]: false,
                    }));
                });
        } else {
            setNameWarnings((prev) => ({
                ...prev,
                [i]: null,
            }));
            setNameValidity((prev) => ({ ...prev, [i]: true }));
        }
    });
}, [substrates.map((s) => s.newSubstrateName).join('|')]);

    return (
        <Box
            sx={{
                flexGrow: 1,
                borderRadius: '11px',
                boxShadow: '0px 4px 10px rgba(100, 84, 31, 0.5)',
                display: 'flex',
                flexDirection: 'column',
                backgroundColor: hasSequenceMatch ? '#e0e0e0' : 'white',
            }}
        >
            {/* Header with collapse toggle */}
            <Box
              sx={{
                backgroundColor: hasSequenceMatch ? '#c0c0c0' : 'secondary.main',
                color: 'black.main',
                padding: '14px 8px',
                display: 'flex',
                justifyContent: 'space-between',
                alignItems: 'center',
                borderTopLeftRadius: '10px',
                borderTopRightRadius: '10px',
              }}
            >
              {/* Left side: domain info + optional warning stacked vertically */}
              <Box sx={{ display: 'flex', flexDirection: 'column' }}>
                <Typography>
                  {`Domain ${parasResult['domain_nr']}: ${domainSynonym} (${parasResult['domain_start']}-${parasResult['domain_end']})`}
                </Typography>

                {hasSequenceMatch && (
                  <Typography
                    variant="body2"
                    sx={{ marginTop: 0.5, fontSize: '0.9em', color: 'black' }}
                  >
                    Sequence already found in dataset ({matchSynonym})
                  </Typography>
                )}
              </Box>

  {/* Right side: expand/collapse button */}
  <IconButton onClick={toggleExpanded} size="small">
    {expanded ? <ExpandLess /> : <ExpandMore />}
  </IconButton>
</Box>

            {/* Collapsible content */}
            <Collapse in={expanded}>
                <Box sx={{padding: 2}}>

                    {/* domain signature */}
                    {parasResult['domain_signature'].length > 0 && (
                        <Box
                            sx={{
                                padding: 1,
                                display: 'flex',
                                justifyContent: 'left',
                            }}
                        >
                            <FaFingerprint size="1.5em" style={{marginRight: '5px'}}/>
                            {parasResult['domain_signature']}
                        </Box>
                    )}

                    {/* spacing */}
                    <Box sx={{height: 16}}/>

                    {/* prediction + substrate selectors + visualizer */}
<Box sx={{ display: 'flex', alignItems: 'flex-start', paddingX: 1, gap: 2, flexDirection: 'column' }}>
    {/* Show known substrates if there's a match */}
    {hasSequenceMatch && (
        <>
            <Typography sx={{ fontWeight: 500 }}>Known substrates</Typography>
            <Box sx={{ display: 'flex', gap: 2, flexWrap: 'wrap' }}>
                {matchedSubstrates.map((substrate, idx) => (
                    <Box key={idx} sx={{ display: 'flex', flexDirection: 'column', alignItems: 'center' }}>
                        <SmileDrawerContainer
                            identifier={`${parasResult['domain_name']}-${parasResult['domain_nr']}-match-${idx}`}
                            smilesStr={substrate.smiles}
                            height={200}
                            width={200}
                        />
                        <Typography variant="caption" sx={{ mt: 1 }}>{substrate.name}</Typography>
                    </Box>
                ))}
            </Box>

            {/* Override checkbox */}
            <FormControlLabel
    control={
        <Checkbox
            checked={overrideSubstrates}
            onChange={(e) => {
                const checked = e.target.checked;
                setOverrideSubstrates(checked);
                if (!checked) {
                    setSubstrates([]);
                    onAnnotationChange?.(domainName, []);
                }
            }}
        />
    }
    label="Make substrate correction"
    sx={{ mt: 2 }}
/>
        </>
    )}

    {/* Render override UI if checked or no match */}
    {(!hasSequenceMatch || overrideSubstrates) && (
        <>


            {/* Render all substrate selector blocks with visualizer */}
            <Box sx={{ display: 'flex', flexDirection: 'column', gap: 2 }}>
                {/* Show PARAS prediction only if no sequence match */}
            {!hasSequenceMatch && selectedPrediction && (
                <Box>
                    <Box
                        sx={{
                            fontWeight: 500,
                            marginBottom: 0.5
                    }}
                        >

                    </Box>

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
                        {`${selectedPrediction['substrate_name']} (${selectedPrediction['probability'].toFixed(2)})`}
                    </Box>
                </Box>
            )}
                <Typography sx={{ fontWeight: 500 }}>Substrate assignments</Typography>
                {substrates.map((sub, i) => (
                    <Box
                        key={i}
                        sx={{
                            display: 'flex',
                            gap: 2,
                            alignItems: 'flex-start',
                            border: '1px solid #ccc',
                            borderRadius: 1,
                            p: 1,
                            minWidth: 500,
                            position: 'relative',
                        }}
                    >
                        {/* Substrate selector box */}
                        <Box sx={{ flex: '1 1 300px' }}>
                            <Box sx={{ fontWeight: 500, marginBottom: 0.5, display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
                                <span>
                                    {sub.substrateName?.trim()
                                        ? `Substrate ${i + 1}: ${sub.substrateName}`
                                        : `Substrate ${i + 1}`}
                                </span>
                                {i > 0 && (
                                    <Button
                                        variant="outlined"
                                        color="error"
                                        size="small"
                                        onClick={() =>
    setSubstrates((prev) => {
        const updated = prev.filter((_, idx) => idx !== i);
        const finalSubs = recalculateAnnotationTypes(updated);
        onAnnotationChange?.(domainName, finalSubs);
        return finalSubs;
    })
}
                                        sx={{ ml: 1, minWidth: 'auto', padding: '2px 8px' }}
                                    >
                                        Remove
                                    </Button>
                                )}
                            </Box>

                            <Autocomplete
                                  options={sortedOptions}
                                  getOptionLabel={(option) => option['name']}
                                  renderInput={(params) => <TextField {...params} variant="outlined" />}
                                  value={sub.selectedSubstrate || null}
                                  onChange={(event, newValue) => updateSubstrateField(i, 'selectedSubstrate', newValue)}
                                  isOptionEqualToValue={(option, value) => option['name'] === value['name']}
                                  clearOnEscape
                                  disabled={sub.useCustomSmiles}
                                  getOptionDisabled={(option) => {
                                    // Disable if this substrate name is already selected in other entries (not counting current index)
                                    return substrates.some((s, idx) =>
                                      idx !== i &&
                                      !s.useCustomSmiles &&
                                      s.selectedSubstrate?.name === option.name
                                    );
                                  }}
                                  PopperComponent={(props) => (
                                      <Popper
                                          {...props}
                                          placement="bottom-start"
                                          modifiers={[{ name: 'flip', enabled: false }]}
                                      />
                                  )}
                                />

                            <FormControlLabel
                                control={
                                    <Checkbox
                                        checked={sub.useCustomSmiles}
                                        onChange={(e) => updateSubstrateField(i, 'useCustomSmiles', e.target.checked)}
                                    />
                                }
                                label="New substrate"
                            />

                            {sub.useCustomSmiles && (
                                <>
                                    <Box sx={{ mt: 1 }}>
                                        <TextField
                                            label="Enter substrate name"
                                            fullWidth
                                            value={sub.newSubstrateName}
                                            onChange={(e) => updateSubstrateField(i, 'newSubstrateName', e.target.value)}
                                            error={Boolean(nameWarnings[i])}
                                            helperText={nameWarnings[i] || ''}
                                        />
                                    </Box>
                                    <Box sx={{ mt: 1 }}>
                                        <SmilesChecker
                                            key={i}
                                            i={i}
                                            customSmiles={sub.customSmiles}
                                            updateSubstrateField={updateSubstrateField}
                                          />
                                    </Box>
                                </>
                            )}
                        </Box>

                        {/* Visualise the substrate */}
                        <Box sx={{ width: 200, height: 200, flexShrink: 0 }}>
                            {/* Prioritize custom SMILES, then selected substrate */}
                            {sub.useCustomSmiles && sub.customSmiles ? (
                                <SmileDrawerContainer
                                    identifier={`${parasResult['domain_name']}-${parasResult['domain_nr']}-custom-${i}`}
                                    smilesStr={sub.customSmiles}
                                    height={200}
                                    width={200}
                                />

                            ) : sub.selectedSubstrate ? (
                                <SmileDrawerContainer
                                    identifier={`${parasResult['domain_name']}-${parasResult['domain_nr']}-sel-${i}`}
                                    smilesStr={sub.selectedSubstrate['smiles']}
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
                                        border: '1px dashed #ccc',
                                    }}
                                >
                                    No substrate selected
                                </Box>
                            )}
                        </Box>
                    </Box>
                ))}

                {/* If the chosen substrates exactly match the given substrates, show message */}

                {showNoUpdateMessage && (
                    <Box
                        sx={{
                            border: '1px solid #aaa',
                            backgroundColor: '#f0f0f0',
                            borderRadius: 1,
                            padding: 2,
                            mb: 2,
                        }}
                    >
                        <Typography sx={{ fontWeight: 500, color: 'green' }}>
                            Selected substrates match known substrates exactly. No update will be made.
                        </Typography>
                    </Box>
                )}

                {/* Add substrate button */}
                <Button
                    variant="outlined"
                    onClick={() =>
                        setSubstrates((prev) => [
                            ...prev,
                            {
                                selectedSubstrate: null,
                                useCustomSmiles: false,
                                customSmiles: '',
                                newSubstrateName: ''
                            },
                        ])
                    }
                    sx={{ alignSelf: 'flex-start', mb: 2 }}
                >
                    Add another substrate
                </Button>
            </Box>
        </>
    )}
</Box>

                </Box>
            </Collapse>
        </Box>
    );
};

export default DomainTile;
