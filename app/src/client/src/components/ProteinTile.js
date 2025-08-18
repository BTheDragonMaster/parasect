import React, {useState, useEffect, useMemo} from 'react';
import {Box, TextField} from '@mui/material';

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
    // Extract first part of protein name (before whitespace)

    const defaultProteinName = useMemo(() => {
        return proteinResult?.protein_name?.split(/\s+/)[0] || "";
    }, [proteinResult]);

    const [proteinName, setProteinName] = useState(defaultProteinName);
    const [proteinExists, setProteinExists] = useState(null); // null = not checked yet

    // 🔍 Check if protein is in dataset
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
        const updatedDomains = {...domainAnnotations};

        if (data && data.length > 0) {
            updatedDomains[domainKey] = {
                name: domainName,
                substrates: data,
                annotationType: annotationType
            };
        } else {
            delete updatedDomains[domainKey];
        }

        setDomainAnnotations(updatedDomains);

        if (onUpdateAnnotation) {
            onUpdateAnnotation(proteinResult["protein_name"], {
                synonym: proteinName,
                sequence: proteinResult["sequence"],
                domains: updatedDomains,
            });
        }
    };


    const handleProteinNameChange = (e) => {
        const newName = e.target.value;
        setProteinName(newName);
        checkProteinInDataset(newName);

        if (onUpdateAnnotation) {
            onUpdateAnnotation(proteinResult["protein_name"], {
                synonym: newName,
                sequence: proteinResult["sequence"],
                domains: domainAnnotations,
            });
        }
    };

    return (
        <Box
            sx={{
                minWidth: '650px',
                maxWidth: '650px',
                borderRadius: '11px',
                boxShadow: '0px 4px 10px rgba(100, 84, 31, 0.5)',
                display: 'flex',
                flexDirection: 'column',
            }}
        >
            {/* header with domain name and location */}
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
                {proteinResult['protein_name']}
            </Box>

            {/* Protein name input field */}
            <Box sx={{px: 2, mb: 1}}>
                <TextField
                    fullWidth
                    label="Protein name"
                    value={proteinName}
                    onChange={handleProteinNameChange}
                    variant="outlined"
                    size="small"
                    error={proteinExists === true}
                    helperText={proteinExists === true ? "This protein name already exists in the dataset." : ""}
                />
            </Box>

            {/* Domains in the protein > 0 */}
            <Box sx={{padding: 2, pt: 1, display: 'flex', flexDirection: 'column', gap: 2}}>
                {proteinResult['results'].map((result, index) => (
                    <DomainTile
                        key={index}
                        domainIndex={index + 1}
                        protein_name={proteinName}
                        result={result}
                        onAnnotationChange={(domainName, data, annotationType) => handleDomainAnnotationChange(index, domainName, data, annotationType)}
                    />

                ))}
            </Box>
        </Box>
    );
};

export default ProteinTile;