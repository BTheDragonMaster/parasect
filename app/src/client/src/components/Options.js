import React, { useState } from "react";
import NumberPicker from "react-widgets/NumberPicker";

const Options = (
    { 
        initVisible,
        selectedSubstrateChoice, 
        setSelectedSubstrateChoice,
        saveActiveSiteSignatures,
        setSaveActiveSiteSignatures,
        saveExtendedSignatures,
        setSaveExtendedSignatures,
        saveAdenylationDomainSequences,
        setSaveAdenylationDomainSequences,
        maxPredictionsToReport,
        setMaxPredictionsToReport,
        numPredictionsToReport,
        setNumPredictionsToReport
    }
) => {
    const [visible, setVisible] = useState(initVisible);

    const handleSubstrateChoiceChange = (event) => {
        setSelectedSubstrateChoice(event.target.value);

        // Update the max number of predictions to report.
        if (event.target.value === "commonSubstrates") {
            setMaxPredictionsToReport(34);
            // setNumPredictionsToReport(34);
        } else {
            setMaxPredictionsToReport(252);
            // setNumPredictionsToReport(252);
        };
    };

    return (
        <div 
            className="control" 
            style={{
                border: "1px solid #dbdbdb", 
                borderRadius: "7.5px", 
                marginBottom: "10px"
            }}
        >
            <div className="panel">

                {/* Panel header. */}
                <div 
                    className="panel-heading"
                    style={{ 
                        display: "flex",
                        justifyContent: "space-between",
                        alignItems: "center"
                    }}
                >
                    <div 
                        className="title is-5" 
                        style={{ margin: "0px" }}
                    >
                        Options
                    </div>
                    <button
                        className="button is-small is-light"
                        onClick={() => setVisible(!visible)}
                    >
                        {visible ? "Hide" : "Show"}
                    </button>
                </div>

                <div>
                    {visible ? (
                        <div>
                            {/* Short description */}
                            <div className="panel-block">
                                <div className="column is-full">
                                    <p>
                                        All selected options for metadata will be included in the results for download.
                                    </p>
                                </div>
                            </div>

                            {/* Checkboxes to select what to save. */}
                            <div 
                                className="panel-block" 
                                style={{ paddingLeft: "20px" }}
                            >
                                <div className="column is-full">
                                    <div className="field has-addons">
                                        <div 
                                            className="control"
                                            style={{ flexDirection: "row", display: "flex", alignItems: "center" }}
                                        >                                            
                                            <div className="checkbox" style={{cursor: "auto"}}>
                                                <input
                                                    id="check1"
                                                    type="checkbox"
                                                    name="check1"
                                                    checked={saveActiveSiteSignatures}
                                                    onChange={() => {setSaveActiveSiteSignatures(!saveActiveSiteSignatures);}}
                                                />
                                                <span style={{marginLeft: "-5px"}}>
                                                    Save active site signatures (Stachelhaus code)
                                                </span> 
                                            </div>
                                        </div>
                                    </div>
                                    <div className="field has-addons">
                                        <div 
                                            className="control"
                                            style={{ flexDirection: "row", display: "flex", alignItems: "center" }}
                                        >   
                                            <div className="checkbox" style={{cursor: "auto"}}>
                                                <input
                                                    id="check2"
                                                    type="checkbox"
                                                    name="check2"
                                                    checked={saveExtendedSignatures}
                                                    onChange={() => {setSaveExtendedSignatures(!saveExtendedSignatures);}}
                                                />
                                                <span style={{ marginLeft: "-5px" }}>
                                                    Save extended signatures (34 amino acid active sites)
                                                </span>
                                            </div>
                                        </div>
                                    </div>
                                    <div className="field has-addons">
                                        <div 
                                            className="control"
                                            style={{ flexDirection: "row", display: "flex", alignItems: "center" }}
                                        >   
                                            <div className="checkbox" style={{cursor: "auto"}}>
                                                <input
                                                    id="check3"
                                                    type="checkbox"
                                                    name="check3"
                                                    checked={saveAdenylationDomainSequences}
                                                    onChange={() => {setSaveAdenylationDomainSequences(!saveAdenylationDomainSequences);}}
                                                />
                                                <span style={{ marginLeft: "-5px" }}>
                                                    Save adenylation domain sequences
                                                </span>
                                            </div>
                                        </div>
                                    </div>
                                </div>
                            </div>

                            {/* Pick model. */}
                            {(selectedSubstrateChoice !== null && setSelectedSubstrateChoice !== null) && (
                                <div 
                                    className="panel-block" 
                                    style={{paddingLeft: "20px"}}
                                >
                                    <div className="column is-full">
                                        <div className="field has-addons">
                                            <div className="control">
                                                <label className="radio">
                                                    <input
                                                        type="radio"
                                                        name="substrateType"
                                                        value="allSubstrates"
                                                        checked={selectedSubstrateChoice === "allSubstrates"}
                                                        onChange={handleSubstrateChoiceChange}
                                                    />
                                                    <span style={{ marginLeft: "5px" }}>
                                                        Use model trained on all substrates
                                                    </span>
                                                </label>
                                            </div>
                                        </div>
                                        <div className="field has-addons">
                                            <div className="control">
                                                <label className="radio">
                                                    <input
                                                        type="radio"
                                                        name="substrateType"
                                                        value="commonSubstrates"
                                                        checked={selectedSubstrateChoice === "commonSubstrates"}
                                                        onChange={handleSubstrateChoiceChange}
                                                    />
                                                    <span style={{ marginLeft: "5px" }}>
                                                        Use model trained on 34 common substrates
                                                    </span>
                                                </label>
                                            </div>
                                        </div>
                                    </div>
                                </div>
                            )}

                            {/* Dropdown for number of predictions to report. */}
                            <div 
                                className="panel-block" 
                                style={{ paddingLeft: "20px", paddingBottom: "20px" }}
                            >
                                <span style={{ marginLeft: "10px", marginRight: "10px" }}>
                                    Number of top predictions to report:
                                    <NumberPicker
                                        value={numPredictionsToReport}
                                        defaultValue={maxPredictionsToReport}
                                        min={1}
                                        max={maxPredictionsToReport}
                                        step={1}  
                                        onChange={setNumPredictionsToReport}
                                    />
                                </span>
                            </div>
                        </div>
                    ) : (
                        <div 
                            className="panel-block" 
                            style={{ padding: "0px", height: "5px" }}
                        />
                )}

                </div>
            </div>
        </div>
    );
};

export default Options;