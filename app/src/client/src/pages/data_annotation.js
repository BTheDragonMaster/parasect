import React, { useState, useEffect, useRef, useMemo } from 'react';
import { toast } from 'react-toastify';
import {
  Box, Button, Divider, TextField, FormControl, FormLabel, RadioGroup,
  FormControlLabel, Radio, CircularProgress, Typography, Input
} from '@mui/material';

const exampleFastaInput = '>dptA\nMDMQSQRLGVTAAQQSVWLAGQLADDHRLYHCAAYLSLTGSIDPRTLGTAVRRTLDETEALRTRFVPQDGELLQILEPGAGQLLLEADFSGDPDPERAAHDWMHAALAAPVRLDRAGTATHALLTLGPSRHLLYFGYHHIALDGYGALLHLRRLAHVYTALSNGDDPGPCPFGPLAGVLTEEAAYRDSDNHRRDGEFWTRSLAGADEAPGLSEREAGALAVPLRRTVELSGERTEKLAASAAATGARWSSLLVAATAAFVRRHAAADDTVIGLPVTARLTGPALRTPCMLANDVPLRLDARLDAPFAALLADTTRAVGTLARHQRFRGEELHRNLGGVGRTAGLARVTVNVLAYVDNIRFGDCRAVVHELSSGPVRDFHINSYGTPGTPDGVQLVFSGNPALYTATDLADHQERFLRFLDAVTADPDLPTGRHRLLSPGTRARLLDDSRGTERPVPRATLPELFAEQARRTPDAPAVQHDGTVLTYRDLHRSVERAAGRLAGLGLRTEDVVALALPKSAESVAILLGIQRAGAAYVPLDPTHPAERLARVLDDTRPRYLVTTGHIDGLSHPTPQLAAADLLREGGPEPAPGRPAPGNAAYIIQTSGSTGRPKGVVVTHEGLATLAADQIRRYRTGPDARVLQFISPGFDVFVSELSMTLLSGGCLVIPPDGLTGRHLADFLAAEAVTTTSLTPGALATMPATDLPHLRTLIVGGEVCPPEIFDQWGRGRDIVNAYGPTETTVEATAWHRDGATHGPVPLGRPTLNRRGYVLDPALEPVPDGTTGELYLAGEGLARGYVAAPGPTAERFVADPFGPPGSRMYRTGDLVRRRSGGMLEFVGRADGQVKLRGFRIELGEVQAALTALPGVRQAGVLIREDRPGDPRLVGYIVPAPGAEPDAGELRAALARTLPPHMVPWALVPLPALPLTSNGKLDRAALPVPAARAGGSGQRPVTPQEKTLCALFADVLGVTEVATDDVFFELGGHSLNGTRLLARIRTEFGTDLTLRDLFAFPTVAGLLPLLDDNGRQHTTPPLPPRPERLPLSHAQQRLWFLDQVEGPSPAYNIPTAVRLEGPLDIPALAVALQDVTNRHEPLRTLLAEDSEGPHQVILPPEAARPELTHSTVAPGDLAAALAEAARRPFDLAGEIPLKAHLFGCGPDDHTLLLLVHHTAGDGASVEVLVRDLAHAYGARRAGDAPHFEPLPLQYADHTLRRRHLLDDPSDSTQLDHWRDALAGLPEQLELPTDHTRPAVPTRRGEAIAFTVPEHTHHTLRAMAQAHGVTVFMVMQAALAALLSRHGAGHDIPLGTPVAGRSDDGTEDLVGFFVNTLVLRNDVSGDPTFAELVSRVRAANLDAYAYQDVPFERLVDVLKPERSLSWHPLFQIMIAYNGPATNDTADGSRFAGLTSRVHAVHTGMSKFDLSFFLTEHADGLGIDGALEFSTDLFTRITAERLVQRYLTVLEQAAGAPDRPISSYELLGDDERALLAQWNDTAHPTPPGTVLDLLESRAARTPDRPAVVENDHVLTYADLHTRANRLARHLITAHGVGPERLVAVALPRSAELLVALLAVLKTGAAYVPLDLTHPAERTAVVLDDCRPAVILTDAGAARELPRRDIPQLRLDEPEVHAAIAEQPGGPVTDRDRTCVTPVSGEHVAYVIYTSGSTGRPKGVAVEHRSLADFVRYSVTAYPGAFDVTLLHSPVTFDLTVTSLFPPLVVGGAIHVADLTEACPPSLAAAGGPTFVKATPSHLPLLTHEATWAASAKVLLVGGEQLLGRELDKWRAGSPEAVVFNDYGPTEATVNCVDFRIDPGQPIGAGPVAIGRPLRNTRVFVLDGGLRAVPVGVVGELHVAGEGLARGYLGQPGLTAERFVACPFGDAGERMYRTGDLVRWRADGMLEFVGRVDDQVKVRGFRIELGEVEAAVAACPGVDRSVVVVREDRPGDRRLVAYVTAAGDEAEGLAPLIVETAAGRLPGYMVPSAVVVLDEIPLTPNGKVDRAALPAPRVAPAAEFRVTGSPREEALCALFAEVLGVERVGVDDGFFDLGGDSILSIQLVARARRAGLEVSVRDVFEHRTVRALAGVVRESGGVAAAVVDSGVGAVERWPVVEWLAERGGGGLGGAVRAFNQSVVVATPAGITWDELRTVLDAVRERHDAWRLRVVDSGDGAWSLRVDAPAPGGEPDWITRHGMASADLEEQVNAVRAAAVEARSRLDPLTGRMVRAVWLDRGPDRRGVLVLVAHHLVVDGVSWRIVLGDLGEAWTQARAGGHVRLDTVGTSLRGWAAALAEQGRHGARATEANLWAQMVHGSDPLVGPRAVDPSVDVFGVVESVGSRASVGVSRALLTEVPSVLGVGVQEVLLAAFGLAVTRWRGRGGSVVVDVEGHGRNEDAVPGADLSRTVGWFTSIYPVRLPLEPAAWDEIRAGGPAVGRTVREIKECLRTLPDQGLGYGILRYLDPENGPALAQHPTPHFGFNYLGRVSVSADAASLDEGDAHADGLGGLVGGRAAADSDEEQWADWVPVSGPFAVGAGQDPVLPVAHAVEFNAITLDTPDGPRLSVTWSWPTTLLSESRIRELARFWDEALEGLVAHARRPDAGGLTPSDLPLVALDHAELEALQADVTGGVHDILPVSPLQEGLLFHSSFAADGVDVYVGQLTFDLTGPVDADHLHAVVESLVTRHDVLRTGYRQAQSGEWIAVVARQVHTPWQYIHTLDTDADTLTNDERWRPFDMTQGPLARFTLARINDTHFRFIVTYHHVILDGWSVAVLIRELFTTYRDTALGRRPEVPYSPPRRDFMAWLAERDQTAAGQAWRSALAGLAEPTVLALGTEGSGVIPEVLEEEISEELTSELVAWARGRGVTVASVVQAAWALVLGRLVGRDDVVFGLTVSGRPAEVAGVEDMVGLFVNTIPLRARMDPAESLGAFVERLQREQTELLEHQHVRLAEVQRWAGHKELFDVGMVFENYPMDSLLQDSLFHGSGLQIDGIQGADATHFALNLAVVPLPAMRFRLGYRPDVFDAGRVRELWGWIVRALECVVCERDVPVSGVDVLGAGERETLLGWGAGAEPGVRALPGAGAGAGAGLVGLFEERVRTDPDAVAVRGAGVEWSYAELNARANAVARWLIGRGVGPERGVGVVMDRGPDVVAMLLAVAKSGGFYLPVDPQWPTERIDWVLADAGIDLAVVGENLAAAVEAVRDCEVVDYAQIARETRLNEQAATDAGDVTDGERVSALLSGHPLYVIYTSGSTGLPKGVVVTHASVGAYLRRGRNAYRGAADGLGHVHSSLAFDLTVTVLFTPLVSGGCVTLGDLDDTANGLGATFLKATPSHLPLLGQLDRVLAPDATLLLGGEALTAGALHHWRTHHPHTTVINAYGPTELTVNCAEYRIPPGHCLPDGPVPIGRPFTGHHLFVLDPALRLTPPDTIGELYVAGDGLARGYLGRPDLTAERFVACPFRSPGERMYRTGDLARWRSDGTLEFIGRADDQVKIRGFRIELGEVEAAVAAHPHVARAIAVVREDRPGDQRLVAYVTGSDPSGLSSAVTDTVAGRLPAYMVPSAVVVLDQIPLTPNGKVDRAALPAPGTASGTTSRAPGTAREEILCTLFADVLGLDQVGVDEDFFDLGGHSLLATRLTSRIRSALGIDLGVRALFKAPTVGRLDQLLQQQTTSLRAPLVARERTGCEPLSFAQQRLWFLHQLEGPNAAYNIPMALRLTGRLDLTALEAALTDVIARHESLRTVIAQDDSGGVWQNILPTDDTRTHLTLDTMPVDAHTLQNRVDEAARHPFDLTTEIPLRATVFRVTDDEHVLLLVLHHIAGDGWSMAPLAHDLSAAYTVRLEHHAPQLPALAVQYADYAAWQRDVLGTENNTSSQLSTQLDYWYSKLEGLPAELTLPTSRVRPAVASHACDRVEFTVPHDVHQGLTALARTQGATVFMVVQAALAALLSRLGAGTDIPIGTPIAGRTDQAMENLIGLFVNTLVLRTDVSGDPTFAELLARVRTTALDAYAHQDIPFERLVEAINPERSLTRHPLFQVMLAFNNTDRRSALDALDAMPGLHARPADVLAVTSPYDLAFSFVETPGSTEMPGILDYATDLFDRSTAEAMTERLVRLLAEIARRPELSVGDIGILSADEVKALSPEAPPAAEELHTSTLPELFEEQVAARGHAVAVVCEGEELSYKELNARANRLARVLMERGAGPERFVGVALPRGLDLIVALLAVTKTGAAYVPLDPEYPTDRLAYMVTDANPTAVVTSTDVHIPLIAPRIELDDEAIRTELAAAPDTAPCVGSGPAHPAYVIYTSGSTGRPKGVVISHANVVRLFTACSDSFDFGPDHVWTLFHSYAFDFSVWEIWGALLHGGRLVVVPFEVTRSPAEFLALLAEQQVTLLSQTPSAFHQLTEAARQEPARCAGLALRHVVFGGEALDPSRLRDWFDLPLGSRPTLVNMYGITETTVHVTVLPLEDRATSLSGSPIGRPLADLQVYVLDERLRPVPPGTVGEMYVAGAGLARGYLGRPALTAERFVADPNSRSGGRLYRTGDLAKVRPDGGLEYVGRGDRQVKIRGFRIELGEIEAALVTHAGVVQAVVLVRDEQTDDQRLVAHVVPALPHRAPTLAELHEHLAATLPAYMVPSAYRTLDELPLTANGKLDRAALAGQWQGGTRTRRLPRTPQEEILCELFADVLRLPAAGADDDFFALGGHSLLATRLLSAVRGTLGVELGIRDLFAAPTPAGLATVLAASGTALPPVTRIDRRPERLPLSFAQRRLWFLSKLEGPSATYNIPVAVRLTGALDVPALRAALGDVTARHESLRTVFPDDGGEPRQLVLPHAEPPFLTHEVTVGEVAEQAASATGYAFDITSDTPLRATLLRVSPEEHVLVVVIHHIAGDGWSMGPLVRDLVTAYRARTRGDAPEYTPLPVQYADYALWQHAVAGDEDAPDGRTARRLGYWREMLAGLPEEHTLPADRPRPVRSSHRGGRVRFELPAGVHRSLLAVARDRRATLFMVVQAALAGLLSRLGAGDDIPIGTPVAGRGDEALDDVVGFFVNTLVLRTNLAGDPSFADLVDRVRTADLDAFAHQDVPFERLVEALAPRRSLARHPLFQIWYTLTNADQDITGQALNALPGLTGDEYPLGASAAKFDLSFTFTEHRTPDGDAAGLSVLLDYSSDLYDHGTAAALGHRLTGFFAALAADPTAPLGTVPLLTDDERDRILGDWGSGTHTPLPPRSVAEQIVRRAALDPDAVAVITAEEELSYRELERLSGETARLLADRGIGRESLVAVALPRTAGLVTTLLGVLRTGAAYLPLDTGYPAERLAHVLSDARPDLVLTHAGLAGRLPAGLAPTVLVDEPQPPAAAAPAVPTSPSGDHLAYVIHTSGSTGRPKGVAIAESSLRAFLADAVRRHDLTPHDRLLAVTTVGFDIAGLELFAPLLAGAAIVLADEDAVRDPASITSLCARHHVTVVQATPSWWRAMLDGAPADAAARLEHVRILVGGEPLPADLARVLTATGAAVTNVYGPTEATIWATAAPLTAGDDRTPGIGTPLDNWRVHILDAALGPVPPGVPGEIHIAGSGLARGYLRRPDLTAERFVANPFAPGERMYRTGDLGRFRPDGTLEHLGRVDDQVKVRGFRIELGDVEAALARHPDVGRAAAAVRPDHRGQGRLVAYVVPRPGTRGPDAGELRETVRELLPDYMVPSAQVTLTTLPHTPNGKLDRAALPAPVFGTPAGRAPATREEKILAGLFADILGLPDVGADSGFFDLGGDSVLSIQLVSRARREGLHITVRDVFEHGTVGALAAAALPAPADDADDTVPGTDVLPSISDDEFEEFELELGLEGEEEQW';

const DataAnnotation = () => {
  const [isLoading, setIsLoading] = useState(false);

  // 'paste' | 'upload' | 'ncbi'
  const [inputMethod, setInputMethod] = useState('paste');

  // 'fasta' | 'gbk' | 'accession'  (accession only valid for 'ncbi')
  const [selectedInputType, setSelectedInputType] = useState('fasta');

  const [selectedInput, setSelectedInput] = useState('');
  const fileInputRef = useRef(null);

  // Keep type consistent with method, and clear stale input on method change
  useEffect(() => {
    if (inputMethod === 'ncbi') {
      setSelectedInputType('accession');
    } else {
      // Leaving NCBI: if type was 'accession', revert to default file text type
      setSelectedInputType(prev => (prev === 'accession' ? 'fasta' : prev));
    }

    // Clear text/file contents when switching methods
    setSelectedInput('');
    if (fileInputRef.current) fileInputRef.current.value = '';
  }, [inputMethod]);

  // Load example
  function handleLoadExample() {
    setSelectedInputType('fasta');
    setSelectedInput(exampleFastaInput);
  }

  // Refresh page
  function handleRefresh() {
    localStorage.removeItem('results');
    window.location.reload();
  }

  // File upload
  const handleFileUpload = async (e) => {
    const file = e.target.files?.[0];
    if (!file) return;

    const reader = new FileReader();
    reader.onload = (event) => {
      setSelectedInput(String(event.target.result || ''));
    };
    reader.readAsText(file);
  };

  // Derive label & accept based on state
  const uploadLabel = useMemo(() => {
    if (inputMethod !== 'upload') return '';
    return `Upload your ${selectedInputType.toUpperCase()} file:`;
  }, [inputMethod, selectedInputType]);

  const fileAccept = useMemo(() => {
    if (selectedInputType === 'gbk') return '.gb,.gbk';
    // default to FASTA extensions
    return '.fa,.fasta,.faa';
  }, [selectedInputType]);

  // Optional: quick validation for NCBI accessions
  const validateAccessions = () => {
    if (inputMethod !== 'ncbi') return true;
    const parts = selectedInput.split(';').map(s => s.trim()).filter(Boolean);
    if (parts.length === 0) {
      toast.error('Please enter at least one accession (separated by ;)');
      return false;
    }
    return true;
  };

  // Submit
  const handleSubmit = async () => {
    if (!selectedInput) {
      toast.error('No input provided');
      return;
    }
    if (!validateAccessions()) return;

    setIsLoading(true);

    // If your backend expects { data: { ... } }, keep this; otherwise send fields directly.
    const payload = {
      selectedInputType,
      selectedInput,
    };

    try {
      const response = await fetch('/api/annotate_data', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ data: payload }), // <- adjust if your API expects raw payload
      });

      if (!response.ok) throw new Error('Network response was not ok!');
      const json = await response.json();

      if (json.status === 'success') {
        const jobId = json.payload?.jobId;
        window.location.href = `/annotation_editor/${jobId}`;
      } else if (json.status === 'warning') {
        toast.warn(json.message);
      } else if (json.status === 'failure') {
        toast.error(json.message);
      }
    } catch (error) {
      console.error('Error:', error);
      toast.error(error.message);
    } finally {
      setIsLoading(false);
    }
  };

  return (
    <>
      <Box
        display='flex'
        flexDirection='column'
        alignItems='flex-start'
        padding={4}
        margin='auto'
      >
        <Typography variant='h4' gutterBottom>
          Upload data for annotation
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
            <FormControlLabel value='paste' control={<Radio />} label='Paste data' />
            <FormControlLabel value='upload' control={<Radio />} label='Upload file' />
            <FormControlLabel value='ncbi' control={<Radio />} label='NCBI accession' />
          </RadioGroup>
        </FormControl>

        {/* input type selection (only for paste/upload) */}
        {(inputMethod === 'paste' || inputMethod === 'upload') && (
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
        )}

        {/* text field or file upload */}
        <Box margin={1} sx={{ width: '100%' }}>
          {inputMethod === 'paste' && (
            <TextField
              label={`Input ${selectedInputType.toUpperCase()} content`}
              multiline
              rows={8}
              fullWidth
              variant='outlined'
              value={selectedInput}
              onChange={(e) => setSelectedInput(e.target.value)}
              margin='normal'
              placeholder={`Paste your ${selectedInputType.toUpperCase()} content here`}
            />
          )}

          {inputMethod === 'ncbi' && (
            <TextField
              label='Input accessions, separated by ;'
              fullWidth
              variant='outlined'
              value={selectedInput}
              onChange={(e) => setSelectedInput(e.target.value)}
              margin='normal'
              placeholder='Paste your NCBI protein accessions here'
            />
          )}

          {inputMethod === 'upload' && (
            <Box width='100%' sx={{ mt: 3, mb: 3 }}>
              <Typography variant='body1' gutterBottom>
                {uploadLabel}
              </Typography>
              <Input
                type='file'
                inputProps={{ accept: fileAccept }}
                onChange={handleFileUpload}
                inputRef={fileInputRef}
              />
            </Box>
          )}

          {/* load example button */}
          {inputMethod === 'paste' && (
            <Button
              variant='text'          // fix: was 'ariant'
              color='primary'
              onClick={handleLoadExample}
              disabled={selectedInputType === 'gbk'}
            >
              Load example input
            </Button>
          )}
        </Box>

        {/* submit and refresh buttons */}
        <Box mt={4} display='flex' justifyContent='flex-start' width='100%' gap={2}>
          <Button variant='contained' color='primary' onClick={handleRefresh}>
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
    </>
  );
};

export default DataAnnotation;