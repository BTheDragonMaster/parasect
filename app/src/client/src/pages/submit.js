import React, { useState } from 'react';
import { toast } from 'react-toastify';
import { Box, Button, TextField, FormControl, FormLabel, RadioGroup, FormControlLabel, Radio, Switch, CircularProgress, Typography, Input, Checkbox, Modal, IconButton } from '@mui/material';
import { MdSettings } from 'react-icons/md';

const exampleFastaInput = '>dptA\nMDMQSQRLGVTAAQQSVWLAGQLADDHRLYHCAAYLSLTGSIDPRTLGTAVRRTLDETEALRTRFVPQDGELLQILEPGAGQLLLEADFSGDPDPERAAHDWMHAALAAPVRLDRAGTATHALLTLGPSRHLLYFGYHHIALDGYGALLHLRRLAHVYTALSNGDDPGPCPFGPLAGVLTEEAAYRDSDNHRRDGEFWTRSLAGADEAPGLSEREAGALAVPLRRTVELSGERTEKLAASAAATGARWSSLLVAATAAFVRRHAAADDTVIGLPVTARLTGPALRTPCMLANDVPLRLDARLDAPFAALLADTTRAVGTLARHQRFRGEELHRNLGGVGRTAGLARVTVNVLAYVDNIRFGDCRAVVHELSSGPVRDFHINSYGTPGTPDGVQLVFSGNPALYTATDLADHQERFLRFLDAVTADPDLPTGRHRLLSPGTRARLLDDSRGTERPVPRATLPELFAEQARRTPDAPAVQHDGTVLTYRDLHRSVERAAGRLAGLGLRTEDVVALALPKSAESVAILLGIQRAGAAYVPLDPTHPAERLARVLDDTRPRYLVTTGHIDGLSHPTPQLAAADLLREGGPEPAPGRPAPGNAAYIIQTSGSTGRPKGVVVTHEGLATLAADQIRRYRTGPDARVLQFISPGFDVFVSELSMTLLSGGCLVIPPDGLTGRHLADFLAAEAVTTTSLTPGALATMPATDLPHLRTLIVGGEVCPPEIFDQWGRGRDIVNAYGPTETTVEATAWHRDGATHGPVPLGRPTLNRRGYVLDPALEPVPDGTTGELYLAGEGLARGYVAAPGPTAERFVADPFGPPGSRMYRTGDLVRRRSGGMLEFVGRADGQVKLRGFRIELGEVQAALTALPGVRQAGVLIREDRPGDPRLVGYIVPAPGAEPDAGELRAALARTLPPHMVPWALVPLPALPLTSNGKLDRAALPVPAARAGGSGQRPVTPQEKTLCALFADVLGVTEVATDDVFFELGGHSLNGTRLLARIRTEFGTDLTLRDLFAFPTVAGLLPLLDDNGRQHTTPPLPPRPERLPLSHAQQRLWFLDQVEGPSPAYNIPTAVRLEGPLDIPALAVALQDVTNRHEPLRTLLAEDSEGPHQVILPPEAARPELTHSTVAPGDLAAALAEAARRPFDLAGEIPLKAHLFGCGPDDHTLLLLVHHTAGDGASVEVLVRDLAHAYGARRAGDAPHFEPLPLQYADHTLRRRHLLDDPSDSTQLDHWRDALAGLPEQLELPTDHTRPAVPTRRGEAIAFTVPEHTHHTLRAMAQAHGVTVFMVMQAALAALLSRHGAGHDIPLGTPVAGRSDDGTEDLVGFFVNTLVLRNDVSGDPTFAELVSRVRAANLDAYAYQDVPFERLVDVLKPERSLSWHPLFQIMIAYNGPATNDTADGSRFAGLTSRVHAVHTGMSKFDLSFFLTEHADGLGIDGALEFSTDLFTRITAERLVQRYLTVLEQAAGAPDRPISSYELLGDDERALLAQWNDTAHPTPPGTVLDLLESRAARTPDRPAVVENDHVLTYADLHTRANRLARHLITAHGVGPERLVAVALPRSAELLVALLAVLKTGAAYVPLDLTHPAERTAVVLDDCRPAVILTDAGAARELPRRDIPQLRLDEPEVHAAIAEQPGGPVTDRDRTCVTPVSGEHVAYVIYTSGSTGRPKGVAVEHRSLADFVRYSVTAYPGAFDVTLLHSPVTFDLTVTSLFPPLVVGGAIHVADLTEACPPSLAAAGGPTFVKATPSHLPLLTHEATWAASAKVLLVGGEQLLGRELDKWRAGSPEAVVFNDYGPTEATVNCVDFRIDPGQPIGAGPVAIGRPLRNTRVFVLDGGLRAVPVGVVGELHVAGEGLARGYLGQPGLTAERFVACPFGDAGERMYRTGDLVRWRADGMLEFVGRVDDQVKVRGFRIELGEVEAAVAACPGVDRSVVVVREDRPGDRRLVAYVTAAGDEAEGLAPLIVETAAGRLPGYMVPSAVVVLDEIPLTPNGKVDRAALPAPRVAPAAEFRVTGSPREEALCALFAEVLGVERVGVDDGFFDLGGDSILSIQLVARARRAGLEVSVRDVFEHRTVRALAGVVRESGGVAAAVVDSGVGAVERWPVVEWLAERGGGGLGGAVRAFNQSVVVATPAGITWDELRTVLDAVRERHDAWRLRVVDSGDGAWSLRVDAPAPGGEPDWITRHGMASADLEEQVNAVRAAAVEARSRLDPLTGRMVRAVWLDRGPDRRGVLVLVAHHLVVDGVSWRIVLGDLGEAWTQARAGGHVRLDTVGTSLRGWAAALAEQGRHGARATEANLWAQMVHGSDPLVGPRAVDPSVDVFGVVESVGSRASVGVSRALLTEVPSVLGVGVQEVLLAAFGLAVTRWRGRGGSVVVDVEGHGRNEDAVPGADLSRTVGWFTSIYPVRLPLEPAAWDEIRAGGPAVGRTVREIKECLRTLPDQGLGYGILRYLDPENGPALAQHPTPHFGFNYLGRVSVSADAASLDEGDAHADGLGGLVGGRAAADSDEEQWADWVPVSGPFAVGAGQDPVLPVAHAVEFNAITLDTPDGPRLSVTWSWPTTLLSESRIRELARFWDEALEGLVAHARRPDAGGLTPSDLPLVALDHAELEALQADVTGGVHDILPVSPLQEGLLFHSSFAADGVDVYVGQLTFDLTGPVDADHLHAVVESLVTRHDVLRTGYRQAQSGEWIAVVARQVHTPWQYIHTLDTDADTLTNDERWRPFDMTQGPLARFTLARINDTHFRFIVTYHHVILDGWSVAVLIRELFTTYRDTALGRRPEVPYSPPRRDFMAWLAERDQTAAGQAWRSALAGLAEPTVLALGTEGSGVIPEVLEEEISEELTSELVAWARGRGVTVASVVQAAWALVLGRLVGRDDVVFGLTVSGRPAEVAGVEDMVGLFVNTIPLRARMDPAESLGAFVERLQREQTELLEHQHVRLAEVQRWAGHKELFDVGMVFENYPMDSLLQDSLFHGSGLQIDGIQGADATHFALNLAVVPLPAMRFRLGYRPDVFDAGRVRELWGWIVRALECVVCERDVPVSGVDVLGAGERETLLGWGAGAEPGVRALPGAGAGAGAGLVGLFEERVRTDPDAVAVRGAGVEWSYAELNARANAVARWLIGRGVGPERGVGVVMDRGPDVVAMLLAVAKSGGFYLPVDPQWPTERIDWVLADAGIDLAVVGENLAAAVEAVRDCEVVDYAQIARETRLNEQAATDAGDVTDGERVSALLSGHPLYVIYTSGSTGLPKGVVVTHASVGAYLRRGRNAYRGAADGLGHVHSSLAFDLTVTVLFTPLVSGGCVTLGDLDDTANGLGATFLKATPSHLPLLGQLDRVLAPDATLLLGGEALTAGALHHWRTHHPHTTVINAYGPTELTVNCAEYRIPPGHCLPDGPVPIGRPFTGHHLFVLDPALRLTPPDTIGELYVAGDGLARGYLGRPDLTAERFVACPFRSPGERMYRTGDLARWRSDGTLEFIGRADDQVKIRGFRIELGEVEAAVAAHPHVARAIAVVREDRPGDQRLVAYVTGSDPSGLSSAVTDTVAGRLPAYMVPSAVVVLDQIPLTPNGKVDRAALPAPGTASGTTSRAPGTAREEILCTLFADVLGLDQVGVDEDFFDLGGHSLLATRLTSRIRSALGIDLGVRALFKAPTVGRLDQLLQQQTTSLRAPLVARERTGCEPLSFAQQRLWFLHQLEGPNAAYNIPMALRLTGRLDLTALEAALTDVIARHESLRTVIAQDDSGGVWQNILPTDDTRTHLTLDTMPVDAHTLQNRVDEAARHPFDLTTEIPLRATVFRVTDDEHVLLLVLHHIAGDGWSMAPLAHDLSAAYTVRLEHHAPQLPALAVQYADYAAWQRDVLGTENNTSSQLSTQLDYWYSKLEGLPAELTLPTSRVRPAVASHACDRVEFTVPHDVHQGLTALARTQGATVFMVVQAALAALLSRLGAGTDIPIGTPIAGRTDQAMENLIGLFVNTLVLRTDVSGDPTFAELLARVRTTALDAYAHQDIPFERLVEAINPERSLTRHPLFQVMLAFNNTDRRSALDALDAMPGLHARPADVLAVTSPYDLAFSFVETPGSTEMPGILDYATDLFDRSTAEAMTERLVRLLAEIARRPELSVGDIGILSADEVKALSPEAPPAAEELHTSTLPELFEEQVAARGHAVAVVCEGEELSYKELNARANRLARVLMERGAGPERFVGVALPRGLDLIVALLAVTKTGAAYVPLDPEYPTDRLAYMVTDANPTAVVTSTDVHIPLIAPRIELDDEAIRTELAAAPDTAPCVGSGPAHPAYVIYTSGSTGRPKGVVISHANVVRLFTACSDSFDFGPDHVWTLFHSYAFDFSVWEIWGALLHGGRLVVVPFEVTRSPAEFLALLAEQQVTLLSQTPSAFHQLTEAARQEPARCAGLALRHVVFGGEALDPSRLRDWFDLPLGSRPTLVNMYGITETTVHVTVLPLEDRATSLSGSPIGRPLADLQVYVLDERLRPVPPGTVGEMYVAGAGLARGYLGRPALTAERFVADPNSRSGGRLYRTGDLAKVRPDGGLEYVGRGDRQVKIRGFRIELGEIEAALVTHAGVVQAVVLVRDEQTDDQRLVAHVVPALPHRAPTLAELHEHLAATLPAYMVPSAYRTLDELPLTANGKLDRAALAGQWQGGTRTRRLPRTPQEEILCELFADVLRLPAAGADDDFFALGGHSLLATRLLSAVRGTLGVELGIRDLFAAPTPAGLATVLAASGTALPPVTRIDRRPERLPLSFAQRRLWFLSKLEGPSATYNIPVAVRLTGALDVPALRAALGDVTARHESLRTVFPDDGGEPRQLVLPHAEPPFLTHEVTVGEVAEQAASATGYAFDITSDTPLRATLLRVSPEEHVLVVVIHHIAGDGWSMGPLVRDLVTAYRARTRGDAPEYTPLPVQYADYALWQHAVAGDEDAPDGRTARRLGYWREMLAGLPEEHTLPADRPRPVRSSHRGGRVRFELPAGVHRSLLAVARDRRATLFMVVQAALAGLLSRLGAGDDIPIGTPVAGRGDEALDDVVGFFVNTLVLRTNLAGDPSFADLVDRVRTADLDAFAHQDVPFERLVEALAPRRSLARHPLFQIWYTLTNADQDITGQALNALPGLTGDEYPLGASAAKFDLSFTFTEHRTPDGDAAGLSVLLDYSSDLYDHGTAAALGHRLTGFFAALAADPTAPLGTVPLLTDDERDRILGDWGSGTHTPLPPRSVAEQIVRRAALDPDAVAVITAEEELSYRELERLSGETARLLADRGIGRESLVAVALPRTAGLVTTLLGVLRTGAAYLPLDTGYPAERLAHVLSDARPDLVLTHAGLAGRLPAGLAPTVLVDEPQPPAAAAPAVPTSPSGDHLAYVIHTSGSTGRPKGVAIAESSLRAFLADAVRRHDLTPHDRLLAVTTVGFDIAGLELFAPLLAGAAIVLADEDAVRDPASITSLCARHHVTVVQATPSWWRAMLDGAPADAAARLEHVRILVGGEPLPADLARVLTATGAAVTNVYGPTEATIWATAAPLTAGDDRTPGIGTPLDNWRVHILDAALGPVPPGVPGEIHIAGSGLARGYLRRPDLTAERFVANPFAPGERMYRTGDLGRFRPDGTLEHLGRVDDQVKVRGFRIELGDVEAALARHPDVGRAAAAVRPDHRGQGRLVAYVVPRPGTRGPDAGELRETVRELLPDYMVPSAQVTLTTLPHTPNGKLDRAALPAPVFGTPAGRAPATREEKILAGLFADILGLPDVGADSGFFDLGGDSVLSIQLVSRARREGLHITVRDVFEHGTVGALAAAALPAPADDADDTVPGTDVLPSISDDEFEEFELELGLEGEEEQW';


const RadioLabel = ({ imageSrc, label }) => (
    <Box display="flex" alignItems="center">
        <Box
            component="img"
            src={imageSrc}
            alt=""
            sx={{ width: 40, height: 40, marginRight: 1 }}
        />
        <Typography variant="body1">{label}</Typography>
    </Box>
);


const Submit = () => {
    // page state
    const [isLoading, setIsLoading] = useState(false);

    // input method and type
    const [inputMethod, setInputMethod] = useState('paste'); // 'paste' or 'upload'
    const [selectedInputType, setSelectedInputType] = useState('fasta'); // fasta or gbk
    const [selectedInput, setSelectedInput] = useState('');

    // options
    const [selectedModel, setSelectedModel] = useState('parasAllSubstrates'); // parasAllSubstrates, parasCommonSubstrates, or parasect
    const [useStructureGuidedAlignment, setUseStructureGuidedAlignment] = useState(false);

    // SMILES file state (for PARASECT model)
    const [smilesFileContent, setSmilesFileContent] = useState(''); // Stores content of uploaded SMILES file
    const [useOnlyUploadedSubstrates, setUseOnlyUploadedSubstrates] = useState(false); // Checkbox state

    // modal state
    const [openSettingsModal, setOpenSettingsModal] = useState(false);

    // open and close modal handlers
    const handleOpenSettingsModal = () => setOpenSettingsModal(true);
    const handleCloseSettingsModal = () => setOpenSettingsModal(false);

    // load example input
    function handleLoadExample() {
        setSelectedInputType('fasta');
        setSelectedInput(exampleFastaInput);
    };

    // refresh the page
    function handleRefresh() {
        // remove results from local storage
        localStorage.removeItem('results');

        // reload the page
        window.location.reload();
    };

    // Handle file upload
    const handleFileUpload = async (e) => {
        const file = e.target.files[0];
        if (file) {
            const reader = new FileReader();
            reader.onload = function (event) {
                const fileContent = event.target.result;
                setSelectedInput(fileContent); // set the file content into selectedInput
            };
            reader.readAsText(file);
        };
    };

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

    // handle form submission
    const handleSubmit = async () => {
        setIsLoading(true);

        const data = {
            selectedInputType: selectedInputType,
            selectedInput: selectedInput,
            selectedModel: selectedModel,
            useStructureGuidedAlignment: useStructureGuidedAlignment,
            smilesFileContent: smilesFileContent,
            useOnlyUploadedSubstrates: useOnlyUploadedSubstrates,
        };

        try {
            const response = await fetch('/api/submit_raw', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ data })
            });

            if (!response.ok) {
                throw new Error('Network response was not ok!');
            };

            const json = await response.json();

            if (json.status === 'success') {
                const jobId = json.payload['jobId'];
                window.location.href = `/results/${jobId}`;
            } else if (json.status === 'warning') {
                toast.warn(json.message);
            } else if (json.status === 'failure') {
                toast.error(json.message);
            };
        } catch (error) {
            console.error('Error:', error);
            toast.error(error.message);
        };

        setIsLoading(false);
    };

    return (
        <>
            <Box 
                display='flex' 
                flexDirection='column' 
                alignItems='left' 
                padding={4} 
                maxWidth={800} 
                margin='auto'
            >
                <Typography variant='h4' gutterBottom>
                    Submit your data
                </Typography>

                {/* input method selection */}
                <FormControl component='fieldset' margin='normal'>
                    <FormLabel component='legend'>Input method</FormLabel>
                    <RadioGroup
                        row
                        value={inputMethod}
                        onChange={(e) => setInputMethod(e.target.value)}
                    >
                        <FormControlLabel value='paste' control={<Radio />} label='Paste data' />
                        <FormControlLabel value='upload' control={<Radio />} label='Upload file' />
                    </RadioGroup>
                </FormControl>

                {/* input type selection */}
                <FormControl component='fieldset' margin='normal'>
                    <FormLabel component='legend'>Input type</FormLabel>
                    <RadioGroup
                        row
                        value={selectedInputType}
                        onChange={(e) => setSelectedInputType(e.target.value)}
                    >
                        <FormControlLabel value='fasta' control={<Radio />} label='FASTA' />
                        <FormControlLabel value='gbk' control={<Radio />} label='GBK' />
                    </RadioGroup>
                </FormControl>

                {/* text field or file upload */}
                {inputMethod === 'paste' ? (
                    <TextField
                        label='Input data'
                        multiline
                        rows={8}
                        fullWidth
                        variant='outlined'
                        value={selectedInput}
                        onChange={(e) => setSelectedInput(e.target.value)}
                        margin='normal'
                        placeholder='Paste your sequence here'
                    />
                ) : (
                    <Box width='100%' sx={{ mt: 3, mb: 3 }}>
                        <Typography variant='body1' gutterBottom>
                            Upload your {selectedInputType.toUpperCase()} file:
                        </Typography>
                        <Input
                            type='file'
                            inputProps={{ accept: '.fasta,.fa,.gbk' }} // accept FASTA or GBK files
                            onChange={handleFileUpload}
                        />
                    </Box>
                )}

                {/* load example button */}
                {inputMethod === 'paste' && (
                    <Button v
                        ariant='text' 
                        color='primary' 
                        onClick={handleLoadExample}
                        disabled={selectedInputType === 'gbk'}
                    >
                        Load example input
                    </Button>
                )}

                {/* model selection */}
                <FormControl component='fieldset' margin='normal'>
                    <FormLabel component='legend'>Select model</FormLabel>
                    <RadioGroup
                        row
                        value={selectedModel}
                        onChange={(e) => setSelectedModel(e.target.value)}
                    >
                        <FormControlLabel 
                            value='parasAllSubstrates' 
                            control={<Radio />} 
                            label={<RadioLabel imageSrc={'/paras.png'} label={'PARAS (all substrates)'} />} 
                        />
                        <FormControlLabel 
                            value='parasCommonSubstrates' 
                            control={<Radio />} 
                            label={<RadioLabel imageSrc={'/paras.png'} label={'PARAS (common substrates)'} />}
                        />
                        <FormControlLabel 
                            value='parasect' 
                            control={<Radio />} 
                            label={<RadioLabel imageSrc={'/parasect.png'} label={'PARASECT'} />}
                        />
                    </RadioGroup>
                </FormControl>

                {/* Submit and Refresh Buttons */}
                <Box mt={4} display='flex' justifyContent='left' width='100%' gap={2}>
                    <Button
                        variant='contained'
                        color='primary'
                        startIcon={<MdSettings size={20} style={{ fill: 'white' }}/>}
                        onClick={handleOpenSettingsModal}
                    >
                        Settings
                    </Button>
                    <Button 
                        variant='contained' 
                        color='primary' 
                        onClick={handleRefresh}
                    >
                        Refresh
                    </Button>
                    <Button
                        variant='contained'
                        color='secondary'
                        onClick={handleSubmit}
                        disabled={isLoading || !selectedInput}
                    >
                        {isLoading ? <CircularProgress size={24} /> : 'Submit'}
                    </Button>
                </Box>
            </Box>

            {/* Settings Modal */}
            <Modal open={openSettingsModal} onClose={handleCloseSettingsModal}>
                <Box
                    width={600}
                    bgcolor='white.main'
                    p={4}
                    mx='auto'
                    my={10}
                    borderRadius={4}
                    boxShadow={3}
                >
                    <Typography variant='h6' gutterBottom sx={{ mb: 2 }}>
                        Settings
                    </Typography>

                    <Box
                        sx={{
                            display: 'flex',
                            flexDirection: 'column',
                            alignItems: 'left',
                            gap: 3
                        }}
                    >
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

                        {/* SMILES file upload for PARASECT */}
                        <Box width='100%' margin='normal'>
                            <Typography 
                                variant='body1' 
                                gutterBottom
                                style={{
                                    color: selectedModel !== 'parasect' ? 'rgba(0, 0, 0, 0.26)' : 'rgba(0, 0, 0, 0.87)'
                                }}
                            >
                                Upload custom list of substrates (CSV format as 'name,SMILES' per line):
                            </Typography>
                            <Input
                                type='file'
                                inputProps={{ accept: '.tsv' }}
                                onChange={handleSmilesFileUpload}
                                disabled={selectedModel !== 'parasect'}
                            />
                        </Box>

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
                    </Box>

                    <Box mt={2}>
                        <Button 
                            variant='contained' 
                            color='primary'
                            onClick={handleCloseSettingsModal}
                            sx={{ mt: 2 }}
                        >
                            Close
                        </Button>
                    </Box>
                </Box>
            </Modal>
        </>
    );
};

export default Submit;