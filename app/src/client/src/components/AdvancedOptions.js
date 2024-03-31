import React, { useState } from "react";
import DropdownList from "react-widgets/DropdownList";

import InfoPopUp from "./InfoPopUp";

const AdvancedOptions = (
    {
        initVisible,
        smilesSrc,
        setSmilesSrc,
        onlyMakePredictionsUploadedSmiles,
        setOnlyMakePredictionsUploadedSmiles,
        useStructureGuidedProfileAlignment,
        setUseStructureGuidedProfileAlignment,
        separators,
        firstSeparator,
        setFirstSeparator,
        secondSeparator,
        setSecondSeparator,
        thirdSeparator,
        setThirdSeparator
    }
) => {
    const [visible, setVisible] = useState(initVisible);

    const handleFileChange = (event) => {
        const file = event.target.files[0];
        const reader = new FileReader();

        reader.onload = (event) => {
            setSmilesSrc(event.target.result);
        };

        reader.readAsText(file);
    };

    console.log((smilesSrc && setSmilesSrc))

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
                            onClick={() => setVisible(!visible)}
                        >
                            {visible ? "Hide" : "Show"}
                        </button>
                    </div>
                </div>

                <div>
                    {visible ? (
                        <div>

                            {/* Only make predictions for submitted SMILES strings. */}
                            {(smilesSrc !== null && setSmilesSrc !== null) && (
                                <div 
                                    className="panel-block" 
                                    style={{paddingLeft: "20px"}}
                                >
                                    <div className="column is-full">
                                        Make predictions for uploaded SMILES strings:

                                        {/* Upload file button. */}
                                        <div className="field is-grouped is-grouped-left">
                                            <div className="control" style={{width: "100%"}}>
                                                <input
                                                    type="file"
                                                    className="button"
                                                    style={{width: "100%"}}
                                                    onChange={handleFileChange}
                                                />
                                            </div>
                                        </div>

                                        {/* Only make predictions for uploaded SMILES strings. */}
                                        {(onlyMakePredictionsUploadedSmiles !== null && setOnlyMakePredictionsUploadedSmiles !== null && smilesSrc.length > 0) && (
                                            <div className="field has-addons">
                                                <div className="control">
                                                    <div 
                                                        className="checkbox" 
                                                        style={{cursor: "auto", flexDirection: "row", display: "flex", alignItems: "center"}}
                                                    >
                                                        <input
                                                            id="check1"
                                                            type="checkbox"
                                                            name="check1"
                                                            checked={onlyMakePredictionsUploadedSmiles}
                                                            onChange={() => setOnlyMakePredictionsUploadedSmiles(!onlyMakePredictionsUploadedSmiles)}
                                                        />
                                                        <span style={{marginLeft: "-5px"}}>
                                                            Only make predictions for uploaded SMILES strings
                                                        </span>
                                                    </div>
                                                </div>
                                            </div>
                                        )}

                                    </div>
                                </div>
                            )}

                            {/* Use structure-guided profile alignment. */}
                            <div 
                                className="panel-block" 
                                style={{paddingLeft: "20px"}}
                            >
                                <div className="column is-full">
                                    <div className="field has-addons">
                                        <div className="control">
                                            <div 
                                                className="checkbox" 
                                                style={{cursor: "auto", flexDirection: "row", display: "flex", alignItems: "center"}}
                                            >
                                                <input
                                                    id="check1"
                                                    type="checkbox"
                                                    name="check1"
                                                    checked={useStructureGuidedProfileAlignment}
                                                    onChange={() => setUseStructureGuidedProfileAlignment(!useStructureGuidedProfileAlignment)}
                                                    disabled={true}
                                                />
                                                <span style={{marginLeft: "-5px"}}>
                                                    Use structure-guided profile alignment instead of pHMM for active site extraction (⚠️ SLOW)
                                                </span>
                                                <InfoPopUp infoText={"Not implemented."} />
                                            </div>
                                        </div>
                                    </div>
                                </div>
                            </div>

                            {/* Custom header format. */}
                            <div 
                                className="panel-block" 
                                style={{paddingLeft: "20px"}}
                            >   
                                <div className="column is-full">
                                    <div style={{marginBottom: "10px"}}>
                                        <span>
                                            Set custom header format (current format:
                                        </span>
                                        <span style={{fontFamily: "monospace", fontWeight: "bold"}}>
                                            {` identifier${firstSeparator}domain${secondSeparator}1${firstSeparator}88${thirdSeparator}478`}
                                        </span>
                                        <span>
                                            ):
                                        </span>
                                    </div>
                                    <div style={{fontFamily: "monospace", flexDirection: "row", display: "flex", alignItems: "center"}}>
                                        identifier
                                        <DropdownList
                                            data={separators}
                                            defaultValue={firstSeparator}
                                            onChange={(value) => setFirstSeparator(value)}
                                            style={{width: "65px"}}
                                        />
                                        domain
                                        <DropdownList
                                            data={separators}
                                            defaultValue={secondSeparator}
                                            onChange={(value) => setSecondSeparator(value)}
                                            style={{width: "65px"}}
                                        />
                                        1
                                        {firstSeparator}
                                        88
                                        <DropdownList
                                            data={separators}
                                            defaultValue={thirdSeparator}
                                            onChange={(value) => setThirdSeparator(value)}
                                            style={{width: "65px"}}
                                        />
                                        478
                                    </div>
                                </div>
                            </div>
                        </div>
                    ) : (
                        <div 
                            className="panel-block" 
                            style={{padding: "0px", height: "5px"}}
                        />
                    )}
                </div>

            </div>
        </div>
    );
};

export default AdvancedOptions;