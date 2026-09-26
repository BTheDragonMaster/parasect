import React, {useState, useEffect, useMemo} from 'react';
import {Box, Checkbox, FormControlLabel, TextField} from '@mui/material';

import DomainTile from '../components/DomainTile';

/**
 * Component to display the results of the prediction.
 *
 * @param {Object} props - The component props.
 * @param {Object} props.proteinResult - The protein result object.
 * @returns {React.ReactElement} - The result tile component.
 */
const ProteinTile = ({proteinResult, onUpdateAnnotation}) => {
    const [domainAnnotations, setDomainAnnotations] = useState({});
    const [excludedDomains, setExcludedDomains] = useState(() => new Set());
    const [proteinExcluded, setProteinExcluded] = useState(false);

    const defaultProteinName = useMemo(() => {
        return proteinResult?.protein_name?.split(/\s+/)[0] || "";
    }, [proteinResult]);

    const [proteinName, setProteinName] = useState(defaultProteinName);
    const [proteinExists, setProteinExists] = useState(null); // null = not checked yet

    // Check if protein is in dataset
    const checkProteinInDataset = async (name) => {
        try {
            const response = await fetch("/api/check_protein_name", {
                method: "POST",
                headers: {
                    "Content-Type": "application/json",
                },
                body: JSON.stringify({protein_name: name}),
            });

            const result = await response.json();
            setProteinExists(result.protein_in_dataset);
        } catch (error) {
            console.error("Error checking protein name:", error);
            setProteinExists(null); // Reset on error
        }
    };

    useEffect(() => {
        if (proteinResult?.protein_name) {
            const initial = proteinResult.protein_name.split(/\s+/)[0];
            setProteinName(initial);
            checkProteinInDataset(initial);
        }
    }, [proteinResult]);


    const handleDomainAnnotationChange = (domainKey, domainName, data, annotationType) => {
        setDomainAnnotations((prev) => {
            const updatedDomains = {...prev};

            if (data && data.length > 0 && annotationType) {
                updatedDomains[domainKey] = {
                    name: domainName,
                    substrates: data,
                    annotationType: annotationType
                };
            } else {
                delete updatedDomains[domainKey];
            }

            return updatedDomains;
        });
    };

    const handleDomainExcludedChange = (domainKey, excluded) => {
        setExcludedDomains((prev) => {
            const updated = new Set(prev);
            if (excluded) {
                updated.add(domainKey);
            } else {
                updated.delete(domainKey);
            }
            return updated;
        });
    };

    useEffect(() => {
        if (!onUpdateAnnotation) return;

        const domains = proteinExcluded
            ? {}
            : Object.fromEntries(
                Object.entries(domainAnnotations).filter(([key]) => !excludedDomains.has(key))
            );

        onUpdateAnnotation(proteinResult["protein_name"], {
            synonym: proteinName,
            sequence: proteinResult["sequence"],
            domains: domains,
        });
    }, [domainAnnotations, excludedDomains, proteinExcluded, proteinName, proteinResult, onUpdateAnnotation]);


    const handleProteinNameChange = (e) => {
        const newName = e.target.value;
        setProteinName(newName);
        checkProteinInDataset(newName);
    };

    return (
        <Box
            sx={{
                minWidth: '650px',
                maxWidth: '650px',
                borderRadius: '11px',
                backgroundColor: 'background.paper',
                border: '1px solid',
                borderColor: 'divider',
                boxShadow: 3,
                display: 'flex',
                flexDirection: 'column',
                overflow: 'hidden',
            }}
        >
            {/* header with domain name and location */}
            <Box
                sx={{
                    backgroundColor: 'accent.main',
                    color: 'accent.contrastText',
                    fontWeight: 600,
                    padding: '21px 16px',
                    display: 'flex',
                    marginBottom: 1,
                    justifyContent: 'space-between',
                    alignItems: 'center',
                    gap: 2,
                }}
            >
                <Box sx={{overflowWrap: 'anywhere'}}>{proteinResult['protein_name']}</Box>
                <FormControlLabel
                    control={
                        <Checkbox
                            size="small"
                            checked={proteinExcluded}
                            onChange={(e) => setProteinExcluded(e.target.checked)}
                            sx={{color: 'inherit', '&.Mui-checked': {color: 'inherit'}}}
                        />
                    }
                    label="Exclude protein"
                    sx={{m: 0, flexShrink: 0, whiteSpace: 'nowrap'}}
                />
            </Box>

            {proteinExcluded && (
                <Box sx={{px: 2, mt: 1, color: 'text.secondary', fontStyle: 'italic'}}>
                    This protein and all its domains are excluded from the submission.
                </Box>
            )}

            {/* Protein name input field */}
            <Box sx={{px: 2, mb: 1, mt: 1, opacity: proteinExcluded ? 0.5 : 1}}>
                <TextField
                    fullWidth
                    label="Protein name"
                    value={proteinName}
                    onChange={handleProteinNameChange}
                    disabled={proteinExcluded}
                    variant="outlined"
                    size="small"
                    error={proteinExists === true}
                    helperText={proteinExists === true ? "This protein name already exists in the dataset." : ""}
                />
            </Box>

            {/* Domains in the protein > 0 */}
            <Box
                sx={{
                    padding: 2, pt: 1, display: 'flex', flexDirection: 'column', gap: 2,
                    opacity: proteinExcluded ? 0.5 : 1,
                    pointerEvents: proteinExcluded ? 'none' : 'auto',
                }}
                aria-disabled={proteinExcluded}
            >
                {proteinResult['results'].map((result, index) => (
                    <DomainTile
                        key={index}
                        domainIndex={index + 1}
                        protein_name={proteinName}
                        result={result}
                        onAnnotationChange={(domainName, data, annotationType) => handleDomainAnnotationChange(String(index), domainName, data, annotationType)}
                        excluded={excludedDomains.has(String(index))}
                        onExcludedChange={(excluded) => handleDomainExcludedChange(String(index), excluded)}
                    />

                ))}
            </Box>
        </Box>
    );
};

export default ProteinTile;