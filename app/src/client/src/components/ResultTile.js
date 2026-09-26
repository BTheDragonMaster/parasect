import React, { useState } from 'react';
import { useNavigate } from 'react-router-dom';
import { toast } from 'react-toastify';
import { Box, Button, MenuItem, Select, FormControl } from '@mui/material';
import { FaFingerprint, FaCopy, FaProjectDiagram } from 'react-icons/fa';

import { SIGNATURE_LENGTH } from './DistanceBadge';

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
    const navigate = useNavigate();
    const extendedSignature = result['domain_extended_signature'];

    // same naming as the reference domains, e.g. "Q04747.3.A2"
    const showInNetwork = () => {
        const params = new URLSearchParams({
            signature: extendedSignature,
            name: `${result['domain_name']}.A${result['domain_nr']}`,
        });
        navigate(`/network?${params}`);
    };
    
    return (
        <Box
            sx={{
                minWidth: '350px',
                maxWidth: '350px',
                borderRadius: '11px',
                backgroundColor: 'background.paper',
                border: '1px solid',
                borderColor: 'divider',
                boxShadow: 2,
                overflow: 'hidden',
            }}
        >

            {/* header with domain name and location */}
            <Box
                sx={{
                    backgroundColor: 'accent.main',
                    color: 'accent.contrastText',
                    padding: '14px 8px',
                    display: 'flex',
                    marginBottom: 1,
                    borderTopLeftRadius: '10px',
                    borderTopRightRadius: '10px',
                    justifyContent: 'center',
                }}
            >
                {`${result['domain_name']} domain ${result['domain_nr']} (${result['domain_start']}-${result['domain_end']})`}
            </Box>

            {/* domain signature, only visualize when length > 0 */}
            {result['domain_signature'].length > 0 &&
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
            }
            
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
                        <FaCopy style={{ marginRight: '5px', fill: 'currentColor' }} />
                        Sequence
                    </Button>
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
                            borderBottom: '1px solid white',
                        }}
                        disabled={result['domain_signature'].length === 0}
                    >
                        <FaCopy style={{ marginRight: '5px', fill: 'currentColor' }} />
                        Signature
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
                            navigator.clipboard.writeText(selectedPrediction['substrate_smiles']);
                            toast.success('Copied the substrate SMILES to clipboard!');
                        }}
                        sx={{ 
                            flexGrow: 1, 
                            width: '50%', 
                            borderRadius: '0',
                            borderBottom: '1px solid white',
                        }}
                        disabled={selectedPrediction['substrate_smiles'].length === 0}
                    >
                        <FaCopy style={{ marginRight: '5px', fill: 'currentColor' }} />
                        SMILES
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
                            borderBottom: '1px solid white',
                        }}
                        disabled={result['domain_extended_signature'].length === 0}
                    >
                        <FaCopy style={{ marginRight: '5px', fill: 'currentColor' }} />
                        Ext. signature
                    </Button>
                </Box>
                <Button
                    variant='contained'
                    color='primary'
                    fullWidth
                    onClick={showInNetwork}
                    sx={{
                        borderRadius: '0',
                        borderBottomLeftRadius: '10px',
                        borderBottomRightRadius: '10px',
                    }}
                    disabled={extendedSignature.length !== SIGNATURE_LENGTH}
                >
                    <FaProjectDiagram style={{ marginRight: '5px', fill: 'currentColor' }} />
                    Show in network
                </Button>
            </Box>

        </Box>
    );
};

export default ResultTile;
