import React from "react";

import Input from "./Input";

const Options = () => {
    return (
        <div className="control" style={{border: "1px solid #dbdbdb", borderRadius: "7.5px", marginBottom: "10px"}}>
            <div className="panel">

                {/* Panel header. */}
                <div className="panel-heading">
                    <div className="title is-5">
                        Options
                    </div>
                </div>

                {/* Checkboxes. */}
                <div className="panel-block" style={{paddingLeft: "20px"}}>
                    <div className="column is-full">
                        <div className="field has-addons">
                            <div className="control">
                                <div className="checkbox">
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
                                <div className="checkbox">
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
                                <div className="checkbox">
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

                {/* Radio buttons. */}
                <div className="panel-block" style={{paddingLeft: "20px"}}>
                </div>

                {/* Dropdown. */}
                <div className="panel-block" style={{paddingLeft: "20px"}}>
                </div>

            </div>
        </div>
    );
};

const AdvancedOptions = () => {
    return (
        <div className="control" style={{border: "1px solid #dbdbdb", borderRadius: "7.5px", marginBottom: "10px"}}>
            <div className="panel">

                {/* Panel header. */}
                <div className="panel-heading">
                    <div className="title is-5">
                        Advanced options
                    </div>
                </div>

                {/* Checkboxes. */}
                <div className="panel-block" style={{paddingLeft: "20px"}}>
                </div>

                {/* More checkboxes. */}
                <div className="panel-block" style={{paddingLeft: "20px"}}>
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