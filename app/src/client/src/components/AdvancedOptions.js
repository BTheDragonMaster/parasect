import React, { useState } from "react";
import DropdownList from "react-widgets/DropdownList";

import InfoPopUp from "./InfoPopUp";

const AdvancedOptions = (
    {
        initVisible,
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
                                                />
                                                <span style={{marginLeft: "-5px"}}>
                                                    Use structure-guided profile alignment instead of pHMM for active site extraction (⚠️ SLOW)
                                                </span>
                                                <InfoPopUp infoText={null} />
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