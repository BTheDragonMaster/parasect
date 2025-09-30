import React, { useState } from 'react';
import { toast } from 'react-toastify';
import { Box, Button, Divider, TextField, FormControl, FormLabel, RadioGroup, FormControlLabel, Radio, CircularProgress, Typography, Input } from '@mui/material';
import { MdSettings, MdBugReport } from 'react-icons/md';

import SettingsModal from '../components/SettingsModal';

const exampleFastaInput = '>dptA\nMDMQSQRLGVTAAQQSVWLAGQLADDHRLYHCAAYLSLTGSIDPRTLGTAVRRTLDETEALRTRFVPQDGELLQILEPGAGQLLLEADFSGDPDPERAAHDWMHAALAAPVRLDRAGTATHALLTLGPSRHLLYFGYHHIALDGYGALLHLRRLAHVYTALSNGDDPGPCPFGPLAGVLTEEAAYRDSDNHRRDGEFWTRSLAGADEAPGLSEREAGALAVPLRRTVELSGERTEKLAASAAATGARWSSLLVAATAAFVRRHAAADDTVIGLPVTARLTGPALRTPCMLANDVPLRLDARLDAPFAALLADTTRAVGTLARHQRFRGEELHRNLGGVGRTAGLARVTVNVLAYVDNIRFGDCRAVVHELSSGPVRDFHINSYGTPGTPDGVQLVFSGNPALYTATDLADHQERFLRFLDAVTADPDLPTGRHRLLSPGTRARLLDDSRGTERPVPRATLPELFAEQARRTPDAPAVQHDGTVLTYRDLHRSVERAAGRLAGLGLRTEDVVALALPKSAESVAILLGIQRAGAAYVPLDPTHPAERLARVLDDTRPRYLVTTGHIDGLSHPTPQLAAADLLREGGPEPAPGRPAPGNAAYIIQTSGSTGRPKGVVVTHEGLATLAADQIRRYRTGPDARVLQFISPGFDVFVSELSMTLLSGGCLVIPPDGLTGRHLADFLAAEAVTTTSLTPGALATMPATDLPHLRTLIVGGEVCPPEIFDQWGRGRDIVNAYGPTETTVEATAWHRDGATHGPVPLGRPTLNRRGYVLDPALEPVPDGTTGELYLAGEGLARGYVAAPGPTAERFVADPFGPPGSRMYRTGDLVRRRSGGMLEFVGRADGQVKLRGFRIELGEVQAALTALPGVRQAGVLIREDRPGDPRLVGYIVPAPGAEPDAGELRAALARTLPPHMVPWALVPLPALPLTSNGKLDRAALPVPAARAGGSGQRPVTPQEKTLCALFADVLGVTEVATDDVFFELGGHSLNGTRLLARIRTEFGTDLTLRDLFAFPTVAGLLPLLDDNGRQHTTPPLPPRPERLPLSHAQQRLWFLDQVEGPSPAYNIPTAVRLEGPLDIPALAVALQDVTNRHEPLRTLLAEDSEGPHQVILPPEAARPELTHSTVAPGDLAAALAEAARRPFDLAGEIPLKAHLFGCGPDDHTLLLLVHHTAGDGASVEVLVRDLAHAYGARRAGDAPHFEPLPLQYADHTLRRRHLLDDPSDSTQLDHWRDALAGLPEQLELPTDHTRPAVPTRRGEAIAFTVPEHTHHTLRAMAQAHGVTVFMVMQAALAALLSRHGAGHDIPLGTPVAGRSDDGTEDLVGFFVNTLVLRNDVSGDPTFAELVSRVRAANLDAYAYQDVPFERLVDVLKPERSLSWHPLFQIMIAYNGPATNDTADGSRFAGLTSRVHAVHTGMSKFDLSFFLTEHADGLGIDGALEFSTDLFTRITAERLVQRYLTVLEQAAGAPDRPISSYELLGDDERALLAQWNDTAHPTPPGTVLDLLESRAARTPDRPAVVENDHVLTYADLHTRANRLARHLITAHGVGPERLVAVALPRSAELLVALLAVLKTGAAYVPLDLTHPAERTAVVLDDCRPAVILTDAGAARELPRRDIPQLRLDEPEVHAAIAEQPGGPVTDRDRTCVTPVSGEHVAYVIYTSGSTGRPKGVAVEHRSLADFVRYSVTAYPGAFDVTLLHSPVTFDLTVTSLFPPLVVGGAIHVADLTEACPPSLAAAGGPTFVKATPSHLPLLTHEATWAASAKVLLVGGEQLLGRELDKWRAGSPEAVVFNDYGPTEATVNCVDFRIDPGQPIGAGPVAIGRPLRNTRVFVLDGGLRAVPVGVVGELHVAGEGLARGYLGQPGLTAERFVACPFGDAGERMYRTGDLVRWRADGMLEFVGRVDDQVKVRGFRIELGEVEAAVAACPGVDRSVVVVREDRPGDRRLVAYVTAAGDEAEGLAPLIVETAAGRLPGYMVPSAVVVLDEIPLTPNGKVDRAALPAPRVAPAAEFRVTGSPREEALCALFAEVLGVERVGVDDGFFDLGGDSILSIQLVARARRAGLEVSVRDVFEHRTVRALAGVVRESGGVAAAVVDSGVGAVERWPVVEWLAERGGGGLGGAVRAFNQSVVVATPAGITWDELRTVLDAVRERHDAWRLRVVDSGDGAWSLRVDAPAPGGEPDWITRHGMASADLEEQVNAVRAAAVEARSRLDPLTGRMVRAVWLDRGPDRRGVLVLVAHHLVVDGVSWRIVLGDLGEAWTQARAGGHVRLDTVGTSLRGWAAALAEQGRHGARATEANLWAQMVHGSDPLVGPRAVDPSVDVFGVVESVGSRASVGVSRALLTEVPSVLGVGVQEVLLAAFGLAVTRWRGRGGSVVVDVEGHGRNEDAVPGADLSRTVGWFTSIYPVRLPLEPAAWDEIRAGGPAVGRTVREIKECLRTLPDQGLGYGILRYLDPENGPALAQHPTPHFGFNYLGRVSVSADAASLDEGDAHADGLGGLVGGRAAADSDEEQWADWVPVSGPFAVGAGQDPVLPVAHAVEFNAITLDTPDGPRLSVTWSWPTTLLSESRIRELARFWDEALEGLVAHARRPDAGGLTPSDLPLVALDHAELEALQADVTGGVHDILPVSPLQEGLLFHSSFAADGVDVYVGQLTFDLTGPVDADHLHAVVESLVTRHDVLRTGYRQAQSGEWIAVVARQVHTPWQYIHTLDTDADTLTNDERWRPFDMTQGPLARFTLARINDTHFRFIVTYHHVILDGWSVAVLIRELFTTYRDTALGRRPEVPYSPPRRDFMAWLAERDQTAAGQAWRSALAGLAEPTVLALGTEGSGVIPEVLEEEISEELTSELVAWARGRGVTVASVVQAAWALVLGRLVGRDDVVFGLTVSGRPAEVAGVEDMVGLFVNTIPLRARMDPAESLGAFVERLQREQTELLEHQHVRLAEVQRWAGHKELFDVGMVFENYPMDSLLQDSLFHGSGLQIDGIQGADATHFALNLAVVPLPAMRFRLGYRPDVFDAGRVRELWGWIVRALECVVCERDVPVSGVDVLGAGERETLLGWGAGAEPGVRALPGAGAGAGAGLVGLFEERVRTDPDAVAVRGAGVEWSYAELNARANAVARWLIGRGVGPERGVGVVMDRGPDVVAMLLAVAKSGGFYLPVDPQWPTERIDWVLADAGIDLAVVGENLAAAVEAVRDCEVVDYAQIARETRLNEQAATDAGDVTDGERVSALLSGHPLYVIYTSGSTGLPKGVVVTHASVGAYLRRGRNAYRGAADGLGHVHSSLAFDLTVTVLFTPLVSGGCVTLGDLDDTANGLGATFLKATPSHLPLLGQLDRVLAPDATLLLGGEALTAGALHHWRTHHPHTTVINAYGPTELTVNCAEYRIPPGHCLPDGPVPIGRPFTGHHLFVLDPALRLTPPDTIGELYVAGDGLARGYLGRPDLTAERFVACPFRSPGERMYRTGDLARWRSDGTLEFIGRADDQVKIRGFRIELGEVEAAVAAHPHVARAIAVVREDRPGDQRLVAYVTGSDPSGLSSAVTDTVAGRLPAYMVPSAVVVLDQIPLTPNGKVDRAALPAPGTASGTTSRAPGTAREEILCTLFADVLGLDQVGVDEDFFDLGGHSLLATRLTSRIRSALGIDLGVRALFKAPTVGRLDQLLQQQTTSLRAPLVARERTGCEPLSFAQQRLWFLHQLEGPNAAYNIPMALRLTGRLDLTALEAALTDVIARHESLRTVIAQDDSGGVWQNILPTDDTRTHLTLDTMPVDAHTLQNRVDEAARHPFDLTTEIPLRATVFRVTDDEHVLLLVLHHIAGDGWSMAPLAHDLSAAYTVRLEHHAPQLPALAVQYADYAAWQRDVLGTENNTSSQLSTQLDYWYSKLEGLPAELTLPTSRVRPAVASHACDRVEFTVPHDVHQGLTALARTQGATVFMVVQAALAALLSRLGAGTDIPIGTPIAGRTDQAMENLIGLFVNTLVLRTDVSGDPTFAELLARVRTTALDAYAHQDIPFERLVEAINPERSLTRHPLFQVMLAFNNTDRRSALDALDAMPGLHARPADVLAVTSPYDLAFSFVETPGSTEMPGILDYATDLFDRSTAEAMTERLVRLLAEIARRPELSVGDIGILSADEVKALSPEAPPAAEELHTSTLPELFEEQVAARGHAVAVVCEGEELSYKELNARANRLARVLMERGAGPERFVGVALPRGLDLIVALLAVTKTGAAYVPLDPEYPTDRLAYMVTDANPTAVVTSTDVHIPLIAPRIELDDEAIRTELAAAPDTAPCVGSGPAHPAYVIYTSGSTGRPKGVVISHANVVRLFTACSDSFDFGPDHVWTLFHSYAFDFSVWEIWGALLHGGRLVVVPFEVTRSPAEFLALLAEQQVTLLSQTPSAFHQLTEAARQEPARCAGLALRHVVFGGEALDPSRLRDWFDLPLGSRPTLVNMYGITETTVHVTVLPLEDRATSLSGSPIGRPLADLQVYVLDERLRPVPPGTVGEMYVAGAGLARGYLGRPALTAERFVADPNSRSGGRLYRTGDLAKVRPDGGLEYVGRGDRQVKIRGFRIELGEIEAALVTHAGVVQAVVLVRDEQTDDQRLVAHVVPALPHRAPTLAELHEHLAATLPAYMVPSAYRTLDELPLTANGKLDRAALAGQWQGGTRTRRLPRTPQEEILCELFADVLRLPAAGADDDFFALGGHSLLATRLLSAVRGTLGVELGIRDLFAAPTPAGLATVLAASGTALPPVTRIDRRPERLPLSFAQRRLWFLSKLEGPSATYNIPVAVRLTGALDVPALRAALGDVTARHESLRTVFPDDGGEPRQLVLPHAEPPFLTHEVTVGEVAEQAASATGYAFDITSDTPLRATLLRVSPEEHVLVVVIHHIAGDGWSMGPLVRDLVTAYRARTRGDAPEYTPLPVQYADYALWQHAVAGDEDAPDGRTARRLGYWREMLAGLPEEHTLPADRPRPVRSSHRGGRVRFELPAGVHRSLLAVARDRRATLFMVVQAALAGLLSRLGAGDDIPIGTPVAGRGDEALDDVVGFFVNTLVLRTNLAGDPSFADLVDRVRTADLDAFAHQDVPFERLVEALAPRRSLARHPLFQIWYTLTNADQDITGQALNALPGLTGDEYPLGASAAKFDLSFTFTEHRTPDGDAAGLSVLLDYSSDLYDHGTAAALGHRLTGFFAALAADPTAPLGTVPLLTDDERDRILGDWGSGTHTPLPPRSVAEQIVRRAALDPDAVAVITAEEELSYRELERLSGETARLLADRGIGRESLVAVALPRTAGLVTTLLGVLRTGAAYLPLDTGYPAERLAHVLSDARPDLVLTHAGLAGRLPAGLAPTVLVDEPQPPAAAAPAVPTSPSGDHLAYVIHTSGSTGRPKGVAIAESSLRAFLADAVRRHDLTPHDRLLAVTTVGFDIAGLELFAPLLAGAAIVLADEDAVRDPASITSLCARHHVTVVQATPSWWRAMLDGAPADAAARLEHVRILVGGEPLPADLARVLTATGAAVTNVYGPTEATIWATAAPLTAGDDRTPGIGTPLDNWRVHILDAALGPVPPGVPGEIHIAGSGLARGYLRRPDLTAERFVANPFAPGERMYRTGDLGRFRPDGTLEHLGRVDDQVKVRGFRIELGDVEAALARHPDVGRAAAAVRPDHRGQGRLVAYVVPRPGTRGPDAGELRETVRELLPDYMVPSAQVTLTTLPHTPNGKLDRAALPAPVFGTPAGRAPATREEKILAGLFADILGLPDVGADSGFFDLGGDSVLSIQLVSRARREGLHITVRDVFEHGTVGALAAAALPAPADDADDTVPGTDVLPSISDDEFEEFELELGLEGEEEQW';

/**
 * Component to submit data to the server.
 * 
 * @param {Object} props - The props of the component.
 * @param {string} props.imageSrc - The path to the image.
 * @param {string} props.label - The label for the radio button.
 * @returns {React.ReactElement} - The submit component.
 */
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


/**
 * Component to submit data to the server.
 * 
 * @returns {React.ReactElement} - The submit component.
 */
const Submit = () => {
    // page state
    const [isLoading, setIsLoading] = useState(false);

    // input method and type
    const [inputMethod, setInputMethod] = useState('upload'); // 'paste' or 'upload'
    const [selectedInputType, setSelectedInputType] = useState('fasta'); // fasta or gbk
    const [selectedInput, setSelectedInput] = useState('');

    // options
    const [selectedModel, setSelectedModel] = useState('parasAllSubstrates'); // parasAllSubstrates, parasCommonSubstrates, or parasect
    const [useStructureGuidedAlignment, setUseStructureGuidedAlignment] = useState(false);

    // SMILES file state (for PARASECT model)
    const [smilesFileContent, setSmilesFileContent] = useState(''); // Stores content of uploaded SMILES file
    const [useOnlyUploadedSubstrates, setUseOnlyUploadedSubstrates] = useState(false); // Checkbox state
    const [uploadedSubstratesFileContentHasHeader, setUploadedSubstratesFileContentHasHeader] = useState(true); // Checkbox state

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

    // handle file upload
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
            uploadedSubstratesFileContentHasHeader: uploadedSubstratesFileContentHasHeader,
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
                margin='auto'
            >
                <Typography variant='h4' gutterBottom>
                    Submit your data
                </Typography>
                <Divider />
                
                {/* input method selection */}
                <FormControl component='fieldset' sx={{ mt: 3 }}>
                    <FormLabel component='legend'>Input method</FormLabel>
                    <RadioGroup
                        row
                        value={inputMethod}
                        onChange={(e) => setInputMethod(e.target.value)}
                    >
                        <FormControlLabel value='upload' control={<Radio />} label='Upload file' />
                        <FormControlLabel value='paste' control={<Radio />} label='Paste data' />
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
                <Box margin={1} >
                    {inputMethod === 'paste' ? (
                        <TextField
                            label='Input adenylation domain data (FASTA or GBK)'
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
                        <Button
                            ariant='text' 
                            color='primary' 
                            onClick={handleLoadExample}
                            disabled={selectedInputType === 'gbk'}
                        >
                            Load example input
                        </Button>
                    )}
                </Box>

                {/* model selection */}
                <FormControl component='fieldset' margin='normal'>
                    <FormLabel component='legend'>Select model</FormLabel>
                    <RadioGroup
                        value={selectedModel}
                        onChange={(e) => setSelectedModel(e.target.value)}
                    >
                        <FormControlLabel 
                            value='parasAllSubstrates' 
                            control={<Radio />} 
                            label={
                                <RadioLabel 
                                    imageSrc={'/paras.png'} 
                                    label={'PARAS (all substrates): predict adenylation domain substrate specificty for a pre-defined list of 223 substrates'}
                                />} 
                        />
                        <FormControlLabel 
                            value='parasCommonSubstrates' 
                            control={<Radio />} 
                            label={
                                <RadioLabel 
                                    imageSrc={'/paras.png'} 
                                    label={'PARAS (common substrates): predict adenylation domain substrate specificty for a pre-defined list of 34 common substrates'}
                                />}
                        />
                        <FormControlLabel 
                            value='parasect' 
                            control={<Radio />} 
                            label={
                                <RadioLabel 
                                    imageSrc={'/parasect.png'} 
                                    label={'PARASECT: predict if a list of pre-defined and/or user-supplied substrates interact with the adenylation domains. This model is traned on fungal and bacterial data.'}
                                />}
                        />
                        <FormControlLabel 
                            value='parasectBacterial' 
                            control={<Radio />} 
                            label={
                                <RadioLabel 
                                    imageSrc={'/parasect.png'} 
                                    label={'PARASECT (bacterial): predict if a list of pre-defined and/or user-supplied substrates interact with the adenylation domains. This model is traned on bacterial data only.'}
                                />
                            }
                        />
                    </RadioGroup>
                </FormControl>

                {/* settings, submit, and refresh buttons */}
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

            {/* settings modal */}
            <SettingsModal
                openSettingsModal={openSettingsModal}
                handleCloseSettingsModal={handleCloseSettingsModal}
                selectedModel={selectedModel}
                useStructureGuidedAlignment={useStructureGuidedAlignment}
                setUseStructureGuidedAlignment={setUseStructureGuidedAlignment}
                smilesFileContent={smilesFileContent}
                setSmilesFileContent={setSmilesFileContent}
                useOnlyUploadedSubstrates={useOnlyUploadedSubstrates}
                setUseOnlyUploadedSubstrates={setUseOnlyUploadedSubstrates}
                uploadedSubstratesFileContentHasHeader={uploadedSubstratesFileContentHasHeader}
                setUploadedSubstratesFileContentHasHeader={setUploadedSubstratesFileContentHasHeader}
            />
        </>
    );
};

export default Submit;
