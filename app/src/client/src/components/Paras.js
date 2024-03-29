import React, { useState } from "react";
import NumberPicker from "react-widgets/NumberPicker";

import Input from "./Input";

const Options = () => {
    const [selectedSubstrateChoice, setSelectedSubstrateChoice] = useState("allSubstrates");
    const [maxPredictionsToReport, setMaxPredictionsToReport] = useState(252);
    const [numPredictionsToReport, setNumPredictionsToReport] = useState(maxPredictionsToReport);

    const handleSubstrateChoiceChange = (event) => {
        setSelectedSubstrateChoice(event.target.value);

        // Update the max number of predictions to report.
        if (event.target.value === "commonSubstrates") {
            setMaxPredictionsToReport(34);
            setNumPredictionsToReport(34);
        } else {
            setMaxPredictionsToReport(252);
            setNumPredictionsToReport(252);
        };
    };

    return (
        <div className="control" style={{border: "1px solid #dbdbdb", borderRadius: "7.5px", marginBottom: "10px"}}>
            <div className="panel">

                {/* Panel header. */}
                <div className="panel-heading">
                    <div className="title is-5">
                        Options
                    </div>
                </div>

                {/* Checkboxes to select what to save. */}
                <div 
                    className="panel-block" 
                    style={{paddingLeft: "20px"}}
                >
                    <div className="column is-full">
                        <div className="field has-addons">
                            <div className="control">
                                <div className="checkbox" style={{cursor: "auto"}}>
                                    <input
                                        id="check1"
                                        type="checkbox"
                                        name="check1"
                                    />
                                    <span style={{marginLeft: "-5px"}}>
                                        Save active site signatures (Stachelhaus code)
                                    </span>
                                </div>
                            </div>
                        </div>
                        <div className="field has-addons">
                            <div className="control">
                                <div className="checkbox" style={{cursor: "auto"}}>
                                    <input
                                        id="check2"
                                        type="checkbox"
                                        name="check2"
                                    />
                                    <span style={{marginLeft: "-5px"}}>
                                        Save extended signatures (34 amino acid active sites)
                                    </span>
                                </div>
                            </div>
                        </div>
                        <div className="field has-addons">
                            <div className="control">
                                <div className="checkbox" style={{cursor: "auto"}}>
                                    <input
                                        id="check3"
                                        type="checkbox"
                                        name="check3"
                                    />
                                    <span style={{marginLeft: "-5px"}}>
                                        Save adenylation domain sequences
                                    </span>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>

                {/* Pick model. */}
                <div 
                    className="panel-block" 
                    style={{paddingLeft: "20px"}}
                >
                    <div className="column is-full">
                        <div className="field has-addons">
                            <div className="control">
                                <label class="radio">
                                    <input
                                        type="radio"
                                        name="substrateType"
                                        value="allSubstrates"
                                        checked={selectedSubstrateChoice === "allSubstrates"}
                                        onChange={handleSubstrateChoiceChange}
                                    />
                                    <span style={{marginLeft: "5px"}}>
                                        Use model trained on all substrates
                                    </span>
                                </label>
                            </div>
                        </div>
                        <div className="field has-addons">
                            <div className="control">
                                <label class="radio">
                                    <input
                                        type="radio"
                                        name="substrateType"
                                        value="commonSubstrates"
                                        checked={selectedSubstrateChoice === "commonSubstrates"}
                                        onChange={handleSubstrateChoiceChange}
                                    />
                                    <span style={{marginLeft: "5px"}}>
                                        Use model trained on 34 common substrates
                                    </span>
                                </label>
                            </div>
                        </div>
                    </div>
                </div>

                {/* Dropdown for number of predictions to report. */}
                <div 
                    className="panel-block" 
                    style={{paddingLeft: "20px", paddingBottom: "20px"}}
                >
                    <span style={{marginLeft: "10px", marginRight: "10px"}}>
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
        </div>
    );
};

const AdvancedOptions = () => {
    const [advancedOptionsVisible, setAdvancedOptionsVisible] = useState(false);

    return (
        <div 
            className="control" 
            style={{border: "1px solid #dbdbdb", borderRadius: "7.5px", marginBottom: "10px"}}
        >
            <div className="panel">

                {/* Panel header. */}
                <div className="panel-heading">
                    <div className="title is-5">
                        Advanced options
                        <button
                            className="button is-small is-light"
                            style={{float: "right", marginLeft: "10px", marginTop: "-5px"}}
                            onClick={() => setAdvancedOptionsVisible(!advancedOptionsVisible)}
                        >
                            {advancedOptionsVisible ? "Hide" : "Show"}
                        </button>
                    </div>
                </div>

                <div>
                    {advancedOptionsVisible ? (
                        <div>
                            <div 
                                className="panel-block" 
                                style={{paddingLeft: "20px"}}
                            >
                                <div className="column is-full">
                                    <div className="field has-addons">
                                        <div className="control">
                                            <div className="checkbox" style={{cursor: "auto"}}>
                                                <input
                                                    id="check1"
                                                    type="checkbox"
                                                    name="check1"
                                                />
                                                <span style={{marginLeft: "-5px"}}>
                                                    Use structure-guided profile alignment instead of pHMM for active site extraction (⚠️ SLOW!)
                                                </span>
                                            </div>
                                        </div>
                                    </div>
                                </div>
                            </div>

                            <div 
                                className="panel-block" 
                                style={{paddingLeft: "20px"}}
                            >
                                <span style={{marginLeft: "10px", marginRight: "10px"}}>
                                    Separators to use in output (default look: identifier|domain_1|88-478):
                                    {/**/}
                                </span>
                            </div>
                        </div>
                    ) : (
                        <div 
                            className="panel-block" 
                            style={{padding: "0px", height: "5px"}}
                        >
                            {/* Content */}
                        </div>
                    )}
                </div>

            </div>
        </div>
    );
};

const Paras = () => {
    return (
        <div className="column is-full">
            <div className="title is-4">
                Paras 
            </div>
            <Input />
            <Options />
            <AdvancedOptions />
        </div>
    );
};

export default Paras;