import React, { useState } from "react";
import { toast } from "react-toastify";

import Input from "./Input";
import Options from "./Options";
import AdvancedOptions from "./AdvancedOptions";
import LoadingOverlay from "./LoadingOverlay";
import Modal from "./Modal";

const Paras = () => {
    // Page state.
    const [isLoading, setIsLoading] = useState(false);
    const [modalActive, setModalActive] = useState(false);

    const toggleModal = () => {
        setModalActive(!modalActive);
    };

    // Input.
    const [selectedInputType, setSelectedInputType] = useState("Fasta"); // Fasta or Genbank
    const [src, setSrc] = useState("");

    // Options.
    const [saveActiveSiteSignatures, setSaveActiveSiteSignatures] = useState(false);
    const [saveExtendedSignatures, setSaveExtendedSignatures] = useState(false);
    const [saveAdenylationDomainSequences, setSaveAdenylationDomainSequences] = useState(false);
    const [selectedSubstrateChoice, setSelectedSubstrateChoice] = useState("allSubstrates"); // allSubstrates or commonSubstrates
    const [maxPredictionsToReport, setMaxPredictionsToReport] = useState(252);
    const [numPredictionsToReport, setNumPredictionsToReport] = useState(maxPredictionsToReport);

    // Advanced options.
    const [useStructureGuidedProfileAlignment, setUseStructureGuidedProfileAlignment] = useState(false);
    const separators = ["|", "_", "-"];
    const [firstSeparator, setFirstSeparator] = useState("|");
    const [secondSeparator, setSecondSeparator] = useState("_");
    const [thirdSeparator, setThirdSeparator] = useState("-");

    function handleLoadExample(setSrc) {
        const example = ">dptA|domain_1\nFAEQARRTPDAPAVQHDGTVLTYRDLHRSVERAAGRLAGLGLRTEDVVALALPKSAESVAILLGIQRAGAAYVPLDPTHPAERLARVLDDTRPRYLVTTGHIDGLSHPTPQLAAADLLREGGPEPAPGRPAPGNAAYIIQTSGSTGRPKGVVVTHEGLATLAADQIRRYRTGPDARVLQFISPGFDVFVSELSMTLLSGGCLVIPPDGLTGRHLADFLAAEAVTTTSLTPGALATMPATDLPHLRTLIVGGEVCPPEIFDQWGRGRDIVNAYGPTETTVEATAWHRDGATHGPVPLGRPTLNRRGYVLDPALEPVPDGTTGELYLAGEGLARGYVAAPGPTAERFVADPFGPPGSRMYRTGDLVRRRSGGMLEFVGRADGQVKLR";
        setSrc(example);
    };

    function handleRefresh() {
        window.location.reload();
    };

    const handleSubmit = async () => {
        setIsLoading(true);

        const data = {
            src: src,
            selectedInputType: selectedInputType,
            saveActiveSiteSignatures: saveActiveSiteSignatures,
            saveExtendedSignatures: saveExtendedSignatures,
            saveAdenylationDomainSequences: saveAdenylationDomainSequences,
            selectedSubstrateChoice: selectedSubstrateChoice,
            maxPredictionsToReport: maxPredictionsToReport,
            numPredictionsToReport: numPredictionsToReport,
            useStructureGuidedProfileAlignment: useStructureGuidedProfileAlignment,
            firstSeparator: firstSeparator,
            secondSeparator: secondSeparator,
            thirdSeparator: thirdSeparator
        };

        try {
            const response = await fetch("/api/submit_paras", {
                method: "POST",
                headers: {"Content-Type": "application/json"},
                body: JSON.stringify({ data })
            });

            if (!response.ok) {
                throw new Error("Network response was not ok!");
            };

            const json = await response.json();

            if (json.status === "success") {
                console.log(json.payload);
                toast.success(json.message);
                toggleModal();
            } else if (json.status === "warning") {
                toast.warn(json.message);
            } else if (json.status === "failure") {
                toast.error(json.message);
            };
        } catch (error) {
            console.error("Error:", error);
            toast.error(error);
        };

        setIsLoading(false);
    };

    return (
        <div className="column is-full">
            {isLoading && (<LoadingOverlay />)}
            <div>
                <Modal 
                    closeModal={toggleModal} 
                    modalState={modalActive} 
                    title="Results"
                >
                    <div>
                        <p>Results will be displayed here.</p>
                    </div>
                </Modal>
                <div className="title is-4">
                    Paras 
                </div>
                <Input 
                    initVisible={true}
                    src={src}
                    setSrc={setSrc}
                    selectedInputType={selectedInputType}
                    setSelectedInputType={setSelectedInputType}
                    handleLoadExample={handleLoadExample}
                    handleRefresh={handleRefresh}
                    handleSubmit={handleSubmit}
                />
                <Options 
                    initVisible={true}
                    selectedSubstrateChoice={selectedSubstrateChoice}
                    setSelectedSubstrateChoice={setSelectedSubstrateChoice}
                    saveActiveSiteSignatures={saveActiveSiteSignatures}
                    setSaveActiveSiteSignatures={setSaveActiveSiteSignatures}
                    saveExtendedSignatures={saveExtendedSignatures}
                    setSaveExtendedSignatures={setSaveExtendedSignatures}
                    saveAdenylationDomainSequences={saveAdenylationDomainSequences}
                    setSaveAdenylationDomainSequences={setSaveAdenylationDomainSequences}
                    maxPredictionsToReport={maxPredictionsToReport}
                    setMaxPredictionsToReport={setMaxPredictionsToReport}
                    numPredictionsToReport={numPredictionsToReport}
                    setNumPredictionsToReport={setNumPredictionsToReport}
                />
                <AdvancedOptions 
                    initVisible={false}
                    useStructureGuidedProfileAlignment={useStructureGuidedProfileAlignment}
                    setUseStructureGuidedProfileAlignment={setUseStructureGuidedProfileAlignment}
                    separators={separators}
                    firstSeparator={firstSeparator}
                    setFirstSeparator={setFirstSeparator}
                    secondSeparator={secondSeparator}
                    setSecondSeparator={setSecondSeparator}
                    thirdSeparator={thirdSeparator}
                    setThirdSeparator={setThirdSeparator}
                />
            </div>
        </div>
    );
};

export default Paras;