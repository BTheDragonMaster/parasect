import React from 'react';
import { Box, Divider, FormControlLabel, Switch, Typography, Input, Checkbox, Modal, IconButton } from '@mui/material';
import { MdClose } from 'react-icons/md';

/**
 * Component to display the settings modal.
 * 
 * @param {Object} props - The props of the component.
 * @param {boolean} props.openSettingsModal - The state of the settings modal.
 * @param {Function} props.handleCloseSettingsModal - The function to close the settings modal.
 * @param {string} props.selectedModel - The selected model.
 * @param {boolean} props.useStructureGuidedAlignment - The state of the structure-guided alignment switch.
 * @param {Function} props.setUseStructureGuidedAlignment - The function to set the state of the structure-guided alignment switch.
 * @param {string} props.smilesFileContent - The content of the SMILES file.
 * @param {Function} props.setSmilesFileContent - The function to set the content of the SMILES file.
 * @param {boolean} props.useOnlyUploadedSubstrates - The state of the use only uploaded substrates checkbox.
 * @param {Function} props.setUseOnlyUploadedSubstrates - The function to set the state of the use only uploaded substrates checkbox.
 * @returns {React.ReactElement} - The settings modal component.
 */
const SettingsModal = ({
    // modal state
    openSettingsModal,
    handleCloseSettingsModal,

    // selected model determines which options are available
    selectedModel,

    // option states
    useStructureGuidedAlignment,
    setUseStructureGuidedAlignment,
    smilesFileContent,
    setSmilesFileContent,
    useOnlyUploadedSubstrates,
    setUseOnlyUploadedSubstrates,
    uploadedSubstratesFileContentHasHeader,
    setUploadedSubstratesFileContentHasHeader,
}) => {

    // handle SMILES file upload for PARASECT
    const handleSmilesFileUpload = async (e) => {
        const file = e.target.files[0];
        if (file) {
            const reader = new FileReader();
            reader.onload = function (event) {
                const fileContent = event.target.result;
                setSmilesFileContent(fileContent); // set the file content into smilesFileContent
            };
            reader.readAsText(file);
        }
    };

    return (
        <Modal open={openSettingsModal} onClose={handleCloseSettingsModal}>
            <Box
                width={800}
                bgcolor='white.main'
                mx='auto'
                my={10}
                borderRadius={4}
                boxShadow={3}
            >
                <Box
                    sx={{
                        backgroundColor: 'secondary.main',
                        borderTopLeftRadius: '14px',
                        borderTopRightRadius: '14px',
                        display: 'flex',
                        justifyContent: 'space-between',
                        alignItems: 'center',
                        padding: 1,
                    }}
                >
                    <Typography 
                        variant='h5' 
                        gutterBottom
                        sx={{ 
                            color: 'black.main', 
                            textAlign: 'center',
                            pl: 2,
                            pt: 2, 
                        }}
                    >
                        Settings
                    </Typography>
                    <IconButton onClick={handleCloseSettingsModal}>
                        <MdClose size={24} />
                    </IconButton>
                </Box>


                <Box
                    sx={{
                        display: 'flex',
                        flexDirection: 'column',
                        alignItems: 'left',
                        gap: 3,
                        p: 4,
                    }}
                >
                    <Box>
                        <Typography variant='body1' gutterBottom>
                            Structure-guided profile alignment uses MUSCLE v3.8.1551 to 
                            align the query sequence to a profile generated from the
                            structure of the query sequence. This can improve the accuracy
                            of the prediction, but may take longer to run. 
                        </Typography>

                        {/* structure-guided alignment switch */}
                        <FormControlLabel
                            control={
                                <Switch
                                    checked={useStructureGuidedAlignment}
                                    onChange={(e) => setUseStructureGuidedAlignment(e.target.checked)}
                                />
                            }
                            label='Use structure-guided profile alignment'
                        />
                    </Box>

                    <Divider />
                    
                    <Box>
                        <Typography variant='body1' gutterBottom>
                            PARASECT predicts if a given substrate has an interaction
                            with the adenylation domain. There is a standard list of 
                            of substrates that PARASECT uses, but you can also upload
                            a custom list of substrates in TSV format.
                        </Typography>

                        {/* SMILES file upload for PARASECT */}
                        <Box width='100%' margin='normal'>
                            <Typography 
                                variant='body1' 
                                gutterBottom
                                style={{
                                    color: (selectedModel !== 'parasect' || selectedModel !== 'parasectBacterial')
                                        ? 'rgba(0, 0, 0, 0.26)'
                                        : 'rgba(0, 0, 0, 0.87)'
                                }}
                            >
                                Upload custom list of substrates (TSV format as 'name\tSMILES' per line):
                            </Typography>
                            <Input
                                type='file'
                                inputProps={{ accept: '.tsv' }}
                                onChange={handleSmilesFileUpload}
                                disabled={!(selectedModel === 'parasect' || selectedModel === 'parasectBacterial')}
                            />
                        </Box>
                        
                        <Box
                            sx={{
                                display: 'flex',
                                flexDirection: 'column',
                                gap: 0,
                                marginTop: 2,
                            }}
                        >
                            {/* option to use only uploaded substrates */}
                            <FormControlLabel
                                control={
                                    <Checkbox
                                        checked={useOnlyUploadedSubstrates}
                                        onChange={(e) => setUseOnlyUploadedSubstrates(e.target.checked)}
                                        disabled={!smilesFileContent}
                                    />
                                }
                                label='Use only uploaded custom substrates'
                            />

                            {/* option to specify if the uploaded substrates file has a header */}
                            <Box
                            // as row
                                sx={{
                                    display: 'flex',
                                    justifyContent: 'left',
                                    alignItems: 'center',
                                    gap: 1,
                                }}
                            >
                                <FormControlLabel
                                    control={
                                        <Checkbox
                                            checked={uploadedSubstratesFileContentHasHeader}
                                            onChange={(e) => setUploadedSubstratesFileContentHasHeader(e.target.checked)}
                                            disabled={!smilesFileContent}
                                        />
                                    }
                                    label='Uploaded custom substrates file has a header line with column names'
                                />
                            </Box>
                        </Box>
                    </Box>
                </Box>
            </Box>
        </Modal>
    );
};

export default SettingsModal;
