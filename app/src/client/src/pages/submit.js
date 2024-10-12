import React, { useState } from "react";
import { toast } from "react-toastify";
import { Box, Button, TextField, FormControl, FormLabel, RadioGroup, FormControlLabel, Radio, Switch, CircularProgress, Typography, Input, Checkbox } from "@mui/material";
import LoadingOverlay from "../components/LoadingOverlay";

const exampleFastaInput = ">dptA\nMDMQSQRLGVTAAQQSVWLAGQLADDHRLYHCAAYLSLTGSIDPRTLGTAVRRTLDETEALRTRFVPQDGELLQILEPGAGQLLLEADFSGDPDPERAAHDWMHAALAAPVRLDRAGTATHALLTLGPSRHLLYFGYHHIALDGYGALLHLRRLAHVYTALSNGDDPGPCPFGPLAGVLTEEAAYRDSDNHRRDGEFWTRSLAGADEAPGLSEREAGALAVPLRRTVELSGERTEKLAASAAATGARWSSLLVAATAAFVRRHAAADDTVIGLPVTARLTGPALRTPCMLANDVPLRLDARLDAPFAALLADTTRAVGTLARHQRFRGEELHRNLGGVGRTAGLARVTVNVLAYVDNIRFGDCRAVVHELSSGPVRDFHINSYGTPGTPDGVQLVFSGNPALYTATDLADHQERFLRFLDAVTADPDLPTGRHRLLSPGTRARLLDDSRGTERPVPRATLPELFAEQARRTPDAPAVQHDGTVLTYRDLHRSVERAAGRLAGLGLRTEDVVALALPKSAESVAILLGIQRAGAAYVPLDPTHPAERLARVLDDTRPRYLVTTGHIDGLSHPTPQLAAADLLREGGPEPAPGRPAPGNAAYIIQTSGSTGRPKGVVVTHEGLATLAADQIRRYRTGPDARVLQFISPGFDVFVSELSMTLLSGGCLVIPPDGLTGRHLADFLAAEAVTTTSLTPGALATMPATDLPHLRTLIVGGEVCPPEIFDQWGRGRDIVNAYGPTETTVEATAWHRDGATHGPVPLGRPTLNRRGYVLDPALEPVPDGTTGELYLAGEGLARGYVAAPGPTAERFVADPFGPPGSRMYRTGDLVRRRSGGMLEFVGRADGQVKLRGFRIELGEVQAALTALPGVRQAGVLIREDRPGDPRLVGYIVPAPGAEPDAGELRAALARTLPPHMVPWALVPLPALPLTSNGKLDRAALPVPAARAGGSGQRPVTPQEKTLCALFADVLGVTEVATDDVFFELGGHSLNGTRLLARIRTEFGTDLTLRDLFAFPTVAGLLPLLDDNGRQHTTPPLPPRPERLPLSHAQQRLWFLDQVEGPSPAYNIPTAVRLEGPLDIPALAVALQDVTNRHEPLRTLLAEDSEGPHQVILPPEAARPELTHSTVAPGDLAAALAEAARRPFDLAGEIPLKAHLFGCGPDDHTLLLLVHHTAGDGASVEVLVRDLAHAYGARRAGDAPHFEPLPLQYADHTLRRRHLLDDPSDSTQLDHWRDALAGLPEQLELPTDHTRPAVPTRRGEAIAFTVPEHTHHTLRAMAQAHGVTVFMVMQAALAALLSRHGAGHDIPLGTPVAGRSDDGTEDLVGFFVNTLVLRNDVSGDPTFAELVSRVRAANLDAYAYQDVPFERLVDVLKPERSLSWHPLFQIMIAYNGPATNDTADGSRFAGLTSRVHAVHTGMSKFDLSFFLTEHADGLGIDGALEFSTDLFTRITAERLVQRYLTVLEQAAGAPDRPISSYELLGDDERALLAQWNDTAHPTPPGTVLDLLESRAARTPDRPAVVENDHVLTYADLHTRANRLARHLITAHGVGPERLVAVALPRSAELLVALLAVLKTGAAYVPLDLTHPAERTAVVLDDCRPAVILTDAGAARELPRRDIPQLRLDEPEVHAAIAEQPGGPVTDRDRTCVTPVSGEHVAYVIYTSGSTGRPKGVAVEHRSLADFVRYSVTAYPGAFDVTLLHSPVTFDLTVTSLFPPLVVGGAIHVADLTEACPPSLAAAGGPTFVKATPSHLPLLTHEATWAASAKVLLVGGEQLLGRELDKWRAGSPEAVVFNDYGPTEATVNCVDFRIDPGQPIGAGPVAIGRPLRNTRVFVLDGGLRAVPVGVVGELHVAGEGLARGYLGQPGLTAERFVACPFGDAGERMYRTGDLVRWRADGMLEFVGRVDDQVKVRGFRIELGEVEAAVAACPGVDRSVVVVREDRPGDRRLVAYVTAAGDEAEGLAPLIVETAAGRLPGYMVPSAVVVLDEIPLTPNGKVDRAALPAPRVAPAAEFRVTGSPREEALCALFAEVLGVERVGVDDGFFDLGGDSILSIQLVARARRAGLEVSVRDVFEHRTVRALAGVVRESGGVAAAVVDSGVGAVERWPVVEWLAERGGGGLGGAVRAFNQSVVVATPAGITWDELRTVLDAVRERHDAWRLRVVDSGDGAWSLRVDAPAPGGEPDWITRHGMASADLEEQVNAVRAAAVEARSRLDPLTGRMVRAVWLDRGPDRRGVLVLVAHHLVVDGVSWRIVLGDLGEAWTQARAGGHVRLDTVGTSLRGWAAALAEQGRHGARATEANLWAQMVHGSDPLVGPRAVDPSVDVFGVVESVGSRASVGVSRALLTEVPSVLGVGVQEVLLAAFGLAVTRWRGRGGSVVVDVEGHGRNEDAVPGADLSRTVGWFTSIYPVRLPLEPAAWDEIRAGGPAVGRTVREIKECLRTLPDQGLGYGILRYLDPENGPALAQHPTPHFGFNYLGRVSVSADAASLDEGDAHADGLGGLVGGRAAADSDEEQWADWVPVSGPFAVGAGQDPVLPVAHAVEFNAITLDTPDGPRLSVTWSWPTTLLSESRIRELARFWDEALEGLVAHARRPDAGGLTPSDLPLVALDHAELEALQADVTGGVHDILPVSPLQEGLLFHSSFAADGVDVYVGQLTFDLTGPVDADHLHAVVESLVTRHDVLRTGYRQAQSGEWIAVVARQVHTPWQYIHTLDTDADTLTNDERWRPFDMTQGPLARFTLARINDTHFRFIVTYHHVILDGWSVAVLIRELFTTYRDTALGRRPEVPYSPPRRDFMAWLAERDQTAAGQAWRSALAGLAEPTVLALGTEGSGVIPEVLEEEISEELTSELVAWARGRGVTVASVVQAAWALVLGRLVGRDDVVFGLTVSGRPAEVAGVEDMVGLFVNTIPLRARMDPAESLGAFVERLQREQTELLEHQHVRLAEVQRWAGHKELFDVGMVFENYPMDSLLQDSLFHGSGLQIDGIQGADATHFALNLAVVPLPAMRFRLGYRPDVFDAGRVRELWGWIVRALECVVCERDVPVSGVDVLGAGERETLLGWGAGAEPGVRALPGAGAGAGAGLVGLFEERVRTDPDAVAVRGAGVEWSYAELNARANAVARWLIGRGVGPERGVGVVMDRGPDVVAMLLAVAKSGGFYLPVDPQWPTERIDWVLADAGIDLAVVGENLAAAVEAVRDCEVVDYAQIARETRLNEQAATDAGDVTDGERVSALLSGHPLYVIYTSGSTGLPKGVVVTHASVGAYLRRGRNAYRGAADGLGHVHSSLAFDLTVTVLFTPLVSGGCVTLGDLDDTANGLGATFLKATPSHLPLLGQLDRVLAPDATLLLGGEALTAGALHHWRTHHPHTTVINAYGPTELTVNCAEYRIPPGHCLPDGPVPIGRPFTGHHLFVLDPALRLTPPDTIGELYVAGDGLARGYLGRPDLTAERFVACPFRSPGERMYRTGDLARWRSDGTLEFIGRADDQVKIRGFRIELGEVEAAVAAHPHVARAIAVVREDRPGDQRLVAYVTGSDPSGLSSAVTDTVAGRLPAYMVPSAVVVLDQIPLTPNGKVDRAALPAPGTASGTTSRAPGTAREEILCTLFADVLGLDQVGVDEDFFDLGGHSLLATRLTSRIRSALGIDLGVRALFKAPTVGRLDQLLQQQTTSLRAPLVARERTGCEPLSFAQQRLWFLHQLEGPNAAYNIPMALRLTGRLDLTALEAALTDVIARHESLRTVIAQDDSGGVWQNILPTDDTRTHLTLDTMPVDAHTLQNRVDEAARHPFDLTTEIPLRATVFRVTDDEHVLLLVLHHIAGDGWSMAPLAHDLSAAYTVRLEHHAPQLPALAVQYADYAAWQRDVLGTENNTSSQLSTQLDYWYSKLEGLPAELTLPTSRVRPAVASHACDRVEFTVPHDVHQGLTALARTQGATVFMVVQAALAALLSRLGAGTDIPIGTPIAGRTDQAMENLIGLFVNTLVLRTDVSGDPTFAELLARVRTTALDAYAHQDIPFERLVEAINPERSLTRHPLFQVMLAFNNTDRRSALDALDAMPGLHARPADVLAVTSPYDLAFSFVETPGSTEMPGILDYATDLFDRSTAEAMTERLVRLLAEIARRPELSVGDIGILSADEVKALSPEAPPAAEELHTSTLPELFEEQVAARGHAVAVVCEGEELSYKELNARANRLARVLMERGAGPERFVGVALPRGLDLIVALLAVTKTGAAYVPLDPEYPTDRLAYMVTDANPTAVVTSTDVHIPLIAPRIELDDEAIRTELAAAPDTAPCVGSGPAHPAYVIYTSGSTGRPKGVVISHANVVRLFTACSDSFDFGPDHVWTLFHSYAFDFSVWEIWGALLHGGRLVVVPFEVTRSPAEFLALLAEQQVTLLSQTPSAFHQLTEAARQEPARCAGLALRHVVFGGEALDPSRLRDWFDLPLGSRPTLVNMYGITETTVHVTVLPLEDRATSLSGSPIGRPLADLQVYVLDERLRPVPPGTVGEMYVAGAGLARGYLGRPALTAERFVADPNSRSGGRLYRTGDLAKVRPDGGLEYVGRGDRQVKIRGFRIELGEIEAALVTHAGVVQAVVLVRDEQTDDQRLVAHVVPALPHRAPTLAELHEHLAATLPAYMVPSAYRTLDELPLTANGKLDRAALAGQWQGGTRTRRLPRTPQEEILCELFADVLRLPAAGADDDFFALGGHSLLATRLLSAVRGTLGVELGIRDLFAAPTPAGLATVLAASGTALPPVTRIDRRPERLPLSFAQRRLWFLSKLEGPSATYNIPVAVRLTGALDVPALRAALGDVTARHESLRTVFPDDGGEPRQLVLPHAEPPFLTHEVTVGEVAEQAASATGYAFDITSDTPLRATLLRVSPEEHVLVVVIHHIAGDGWSMGPLVRDLVTAYRARTRGDAPEYTPLPVQYADYALWQHAVAGDEDAPDGRTARRLGYWREMLAGLPEEHTLPADRPRPVRSSHRGGRVRFELPAGVHRSLLAVARDRRATLFMVVQAALAGLLSRLGAGDDIPIGTPVAGRGDEALDDVVGFFVNTLVLRTNLAGDPSFADLVDRVRTADLDAFAHQDVPFERLVEALAPRRSLARHPLFQIWYTLTNADQDITGQALNALPGLTGDEYPLGASAAKFDLSFTFTEHRTPDGDAAGLSVLLDYSSDLYDHGTAAALGHRLTGFFAALAADPTAPLGTVPLLTDDERDRILGDWGSGTHTPLPPRSVAEQIVRRAALDPDAVAVITAEEELSYRELERLSGETARLLADRGIGRESLVAVALPRTAGLVTTLLGVLRTGAAYLPLDTGYPAERLAHVLSDARPDLVLTHAGLAGRLPAGLAPTVLVDEPQPPAAAAPAVPTSPSGDHLAYVIHTSGSTGRPKGVAIAESSLRAFLADAVRRHDLTPHDRLLAVTTVGFDIAGLELFAPLLAGAAIVLADEDAVRDPASITSLCARHHVTVVQATPSWWRAMLDGAPADAAARLEHVRILVGGEPLPADLARVLTATGAAVTNVYGPTEATIWATAAPLTAGDDRTPGIGTPLDNWRVHILDAALGPVPPGVPGEIHIAGSGLARGYLRRPDLTAERFVANPFAPGERMYRTGDLGRFRPDGTLEHLGRVDDQVKVRGFRIELGDVEAALARHPDVGRAAAAVRPDHRGQGRLVAYVVPRPGTRGPDAGELRETVRELLPDYMVPSAQVTLTTLPHTPNGKLDRAALPAPVFGTPAGRAPATREEKILAGLFADILGLPDVGADSGFFDLGGDSVLSIQLVSRARREGLHITVRDVFEHGTVGALAAAALPAPADDADDTVPGTDVLPSISDDEFEEFELELGLEGEEEQW";


const Submit = () => {
    // page state
    const [isLoading, setIsLoading] = useState(false);

    // input method and type
    const [inputMethod, setInputMethod] = useState("paste"); // 'paste' or 'upload'
    const [selectedInputType, setSelectedInputType] = useState("fasta"); // fasta or gbk
    const [selectedInput, setSelectedInput] = useState("");

    // options
    const [selectedModel, setSelectedModel] = useState("parasAllSubstrates"); // parasAllSubstrates, parasCommonSubstrates, or parasect
    const [useStructureGuidedAlignment, setUseStructureGuidedAlignment] = useState(false);

    // SMILES file state (for PARASECT model)
    const [smilesFileContent, setSmilesFileContent] = useState(""); // Stores content of uploaded SMILES file
    const [useOnlyUploadedSubstrates, setUseOnlyUploadedSubstrates] = useState(false); // Checkbox state

    // load example input
    function handleLoadExample() {
        setSelectedInputType("fasta");
        setSelectedInput(exampleFastaInput);
    };

    // refresh the page
    function handleRefresh() {
        // remove results from local storage
        localStorage.removeItem("results");

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
        }
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
            const response = await fetch("/api/submit_raw", {
                method: "POST",
                headers: { "Content-Type": "application/json" },
                body: JSON.stringify({ data })
            });

            if (!response.ok) {
                throw new Error("Network response was not ok!");
            }

            const json = await response.json();

            if (json.status === "success") {
                localStorage.setItem("results", JSON.stringify(json.payload["results"]));
                toast.success(json.message);
                window.location.href = "/results";
            } else if (json.status === "warning") {
                toast.warn(json.message);
            } else if (json.status === "failure") {
                toast.error(json.message);
            }
        } catch (error) {
            console.error("Error:", error);
            toast.error(error.message);
        };

        setIsLoading(false);
    };

    return (
        <>
            {isLoading && (
                <LoadingOverlay 
                    frame1="paras_loading_1.png"
                    frame2="paras_loading_2.png"
                />
            )}
            <Box display="flex" flexDirection="column" alignItems="center" padding={4} maxWidth={800} margin="auto">
                <Typography variant="h4" gutterBottom>
                    Submit your data
                </Typography>

                {/* input method selection */}
                <FormControl component="fieldset" margin="normal">
                    <FormLabel component="legend">Input method</FormLabel>
                    <RadioGroup
                        row
                        value={inputMethod}
                        onChange={(e) => setInputMethod(e.target.value)}
                    >
                        <FormControlLabel value="paste" control={<Radio />} label="Paste sequence" />
                        <FormControlLabel value="upload" control={<Radio />} label="Upload file" />
                    </RadioGroup>
                </FormControl>

                {/* input type selection */}
                <FormControl component="fieldset" margin="normal">
                    <FormLabel component="legend">Input type</FormLabel>
                    <RadioGroup
                        row
                        value={selectedInputType}
                        onChange={(e) => setSelectedInputType(e.target.value)}
                    >
                        <FormControlLabel value="fasta" control={<Radio />} label="FASTA" />
                        <FormControlLabel value="gbk" control={<Radio />} label="GBK" />
                    </RadioGroup>
                </FormControl>

                {/* text field or file upload */}
                {inputMethod === "paste" ? (
                    <TextField
                        label="Input Data"
                        multiline
                        rows={8}
                        fullWidth
                        variant="outlined"
                        value={selectedInput}
                        onChange={(e) => setSelectedInput(e.target.value)}
                        margin="normal"
                        placeholder="Paste your sequence here"
                    />
                ) : (
                    <Box width="100%" margin="normal">
                        <Typography variant="body1" gutterBottom>
                            Upload your {selectedInputType.toUpperCase()} file:
                        </Typography>
                        <Input
                            type="file"
                            inputProps={{ accept: ".fasta,.fa,.gbk" }} // accept FASTA or GBK files
                            onChange={handleFileUpload}
                        />
                    </Box>
                )}

                {/* load example button */}
                {inputMethod === "paste" && (
                    <Button variant="text" color="primary" onClick={handleLoadExample}>
                        Load example input
                    </Button>
                )}

                {/* model selection */}
                <FormControl component="fieldset" margin="normal">
                    <FormLabel component="legend">Select model</FormLabel>
                    <RadioGroup
                        row
                        value={selectedModel}
                        onChange={(e) => setSelectedModel(e.target.value)}
                    >
                        <FormControlLabel value="parasAllSubstrates" control={<Radio />} label="PARAS All Substrates" />
                        <FormControlLabel value="parasCommonSubstrates" control={<Radio />} label="PARAS Common Substrates" />
                        <FormControlLabel value="parasect" control={<Radio />} label="PARASECT" />
                    </RadioGroup>
                </FormControl>

                {/* structure-guided alignment switch */}
                <FormControlLabel
                    control={
                        <Switch
                            checked={useStructureGuidedAlignment}
                            onChange={(e) => setUseStructureGuidedAlignment(e.target.checked)}
                        />
                    }
                    label="Use structure-guided profile alignment (watch out: slow!)"
                />

                {/* SMILES file upload for PARASECT */}
                {selectedModel === "parasect" && (
                    <Box width="100%" margin="normal">
                        <Typography variant="body1" gutterBottom>
                            Upload your list of substrates (SMILES format, one per line):
                        </Typography>
                        <Input
                            type="file"
                            inputProps={{ accept: ".txt" }} // accept only text files for SMILES
                            onChange={handleSmilesFileUpload}
                        />

                        {/* option to use only uploaded substrates */}
                        <FormControlLabel
                            control={
                                <Checkbox
                                    checked={useOnlyUploadedSubstrates}
                                    onChange={(e) => setUseOnlyUploadedSubstrates(e.target.checked)}
                                    disabled={!smilesFileContent}
                                />
                            }
                            label="Use only uploaded substrates"
                        />
                    </Box>
                )}

                {/* Submit and Refresh Buttons */}
                <Box mt={4} display="flex" justifyContent="space-between" width="100%">
                    <Button variant="outlined" color="secondary" onClick={handleRefresh}>
                        Refresh
                    </Button>
                    <Button
                        variant="contained"
                        color="primary"
                        onClick={handleSubmit}
                        disabled={isLoading || !selectedInput}
                    >
                        {isLoading ? <CircularProgress size={24} /> : "Submit"}
                    </Button>
                </Box>
            </Box>
        </>
    );
};

export default Submit;