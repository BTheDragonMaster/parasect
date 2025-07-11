import React, { useState, useEffect } from 'react';
import { toast } from 'react-toastify';
import { Box, Button, MenuItem, Select, FormControl, TextField } from '@mui/material';
import { FaFingerprint, FaCopy } from 'react-icons/fa';

import SmileDrawerContainer from './SmilesDrawer';
import DomainTile from '../components/DomainTile';

/**
 * Component to display the results of the prediction.
 *
 * @param {Object} props - The component props.
 * @param {Object} props.proteinResult - The protein result object.
 * @returns {React.ReactElement} - The result tile component.
 */
const ProteinTile = ({ proteinResult, onUpdateAnnotation }) => {
    const [domainAnnotations, setDomainAnnotations] = useState({});

    // Extract first part of protein name (before whitespace)
    const initialProteinName = proteinResult['protein_name'].split(/\s+/)[0];

    const [proteinName, setProteinName] = useState(initialProteinName);

    const handleDomainAnnotationChange = (domainKey, data) => {
    const updatedDomains = { ...domainAnnotations };

    if (data && data.length > 0) {
        updatedDomains[domainKey] = data;
    } else {
        delete updatedDomains[domainKey];
    }

    setDomainAnnotations(updatedDomains);

    if (onUpdateAnnotation) {
        onUpdateAnnotation(proteinResult["protein_name"], {
            synonym: proteinName,
            domains: updatedDomains,
        });
    }
};

const handleProteinNameChange = (e) => {
    const newName = e.target.value;
    setProteinName(newName);

    if (onUpdateAnnotation) {
        onUpdateAnnotation(proteinResult["protein_name"], {
            synonym: newName,
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
            <Box sx={{ px: 2, mb: 1 }}>
                <TextField
                    fullWidth
                    label="Protein name"
                    value={proteinName}
                    onChange={handleProteinNameChange}
                    variant="outlined"
                    size="small"
                />
            </Box>

            {/* Domains in the protein > 0 */}
            <Box sx={{ padding: 2, pt: 1, display: 'flex', flexDirection: 'column', gap: 2 }}>
                {proteinResult['results'].map((result, index) => (
                    <DomainTile
                        key={index}
                        result={result}
                        onAnnotationChange={(data) => handleDomainAnnotationChange(index, data)}
                    />
                ))}
            </Box>
        </Box>
    );
};

export default ProteinTile;