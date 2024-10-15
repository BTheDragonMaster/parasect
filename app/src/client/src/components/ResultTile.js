import React, { useState } from 'react';
import { toast } from 'react-toastify';
import { Box, Button, MenuItem, Select, FormControl } from '@mui/material';
import { FaFingerprint, FaCopy } from 'react-icons/fa';

import SmileDrawerContainer from './SmilesDrawer';

/**
 * Component to display the results of the prediction.
 * 
 * @param {Object} props - The component props.
 * @param {Object} props.result - The result object.
 * @returns {React.ReactElement} - The result tile component.
 */
const ResultTile = ({ result }) => {
    const [selectedPrediction, setSelectedPrediction] = useState(result['predictions'][0]);
    
    return (
        <Box
            sx={{
                minWidth: '350px',
                maxWidth: '350px',
                borderRadius: '11px',
                boxShadow: '0px 4px 10px rgba(100, 84, 31, 0.5)',
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
                {`${result['domain_name']} domain ${result['domain_nr']} (${result['domain_start']}-${result['domain_start']})`}
            </Box>

            {/* domain signature */}
            <Box 
                sx={{ 
                    padding: 1,
                    display: 'flex',
                    justifyContent: 'center',
                }}
            >
                <FaFingerprint 
                    size='1.5em' 
                    style={{ marginRight: '5px' }}
                />
                {result['domain_signature']}
            </Box>
            
            {/* drop down for picking substrate prediction, sorted on prediction value */}
            <Box sx={{ padding: 1 }}>
                <FormControl fullWidth>
                    <Select
                        labelId="substrate-select"
                        id="substrate-select"
                        value={selectedPrediction['substrate_name']}
                        onChange={(event) => {
                            const selectedPrediction = result['predictions'].find(prediction => prediction['substrate_name'] === event.target.value);
                            setSelectedPrediction(selectedPrediction);
                        }}
                        sx={{ borderRadius: '0' }}
                        MenuProps={{
                            PaperProps: {
                                style: {
                                    maxHeight: 300,
                                },
                            },
                        }}
                    >
                        {result['predictions'].map((substrate, index) => (
                            <MenuItem 
                                key={index} 
                                value={substrate['substrate_name']}
                            >
                                {index + 1}. {substrate['substrate_name']} ({substrate['probability']})
                            </MenuItem>
                        ))}
                    </Select>
                </FormControl>
            </Box>

            {/* visualize selected substrate */}
            <Box
                sx={{
                    display: 'flex',
                    justifyContent: 'center',
                    alignItems: 'center',
                    padding: 1,
                }}
            >
                <SmileDrawerContainer  
                    identifier={`${result['domain_name']}-${result['domain_nr']}`}
                    smilesStr={selectedPrediction['substrate_smiles']} 
                    height={200}
                    width={200}
                />
            </Box>
            
            { /* copy domain information to clipboard */}
            <Box>
                <Box
                    sx={{
                        display: 'flex',
                        justifyContent: 'center',
                        alignItems: 'center',
                        gap: '1px'
                    }}
                >
                    <Button
                        variant='contained'
                        color='primary'
                        onClick={() => {
                            navigator.clipboard.writeText(selectedPrediction['substrate_smiles']);
                            toast.success('Copied the substrate SMILES to clipboard!');
                        }}
                        sx={{ 
                            flexGrow: 1, 
                            width: '50%', 
                            borderRadius: '0',
                            borderBottom: '1px solid',
                        }}
                    >
                        <FaCopy style={{ marginRight: '5px', fill: 'white' }} />
                        SMILES
                    </Button>
                    <Button
                        variant='contained'
                        color='primary'
                        onClick={() => {
                            navigator.clipboard.writeText(result['domain_sequence']);
                            toast.success('Copied the domain amino acid sequence to clipboard!');
                        }}
                        sx={{ 
                            flexGrow: 1, 
                            width: '50%',
                            borderRadius: '0',
                            borderBottom: '1px solid',
                        }}
                    >
                        <FaCopy style={{ marginRight: '5px', fill: 'white' }} />
                        Sequence
                    </Button>
                </Box>
                <Box
                    sx={{
                        display: 'flex',
                        justifyContent: 'center',
                        alignItems: 'center',
                        gap: '1px'
                    }}
                >
                    <Button
                        variant='contained'
                        color='primary'
                        onClick={() => {
                            navigator.clipboard.writeText(result['domain_signature']);
                            toast.success('Copied the domain signature to clipboard!');
                        }}
                        sx={{ 
                            flexGrow: 1, 
                            width: '50%', 
                            borderRadius: '0',
                            borderBottomLeftRadius: '10px',
                        }}
                    >
                        <FaCopy style={{ marginRight: '5px', fill: 'white' }} />
                        Signature
                    </Button>
                    <Button
                        variant='contained'
                        color='primary'
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
                    >
                        <FaCopy style={{ marginRight: '5px', fill: 'white' }} />
                        Ext. signature
                    </Button>
                </Box>
            </Box>

        </Box>
    );
};

export default ResultTile;
