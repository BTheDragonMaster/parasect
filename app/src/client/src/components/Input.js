import React, { useState } from "react";
import { toast } from "react-toastify";

const FastaInputContainer = ({ src, setSrc }) => {

    const handleFileChange = (event) => {
        const file = event.target.files[0];
        const reader = new FileReader();

        reader.onload = (event) => {
            setSrc(event.target.result);
        };

        reader.readAsText(file);
    };

    return (
        <div className="panel-block">
            <div className="column is-full">

                {/* Text input. */}
                <div className="field" style={{width: "100%"}}>
                    <div className="control">
                        <textarea
                            className="textarea"
                            placeholder="Enter Fasta sequence here or upload Fasta file..."
                            value={src}
                            onChange={(event) => setSrc(event.target.value)}
                        >
                        </textarea>
                    </div>
                </div>

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

            </div>
        </div>
    );
};

const GenbankInputContainer = ({ src, setSrc }) => {

    const handleFileChange = (event) => {
        const file = event.target.files[0];
        const reader = new FileReader();

        reader.onload = (event) => {
            setSrc(event.target.result);
        };

        reader.readAsText(file);
    };

    return (
        <div className="panel-block">
            <div className="column is-full">

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

            </div>
        </div>
    );
};

const Input = (
    { 
        initVisible,
        src, 
        setSrc, 
        selectedInputType, 
        setSelectedInputType,
        handleLoadExample,
        handleRefresh,
        handleSubmit
    }
) => {
    const [visible, setVisible] = useState(initVisible);

    const handleInputTypeChange = (event) => {
        setSelectedInputType(event.target.value);
    };

    return (
        <div className="control" style={{border: "1px solid #dbdbdb", borderRadius: "7.5px", marginBottom: "10px"}}>
            <div className="panel">

                {/* Panel header. */}
                <div className="panel-heading">
                    <div className="title is-5">
                        Input
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
                            {/* Add a radio button to select the input type: Fasta or Genbank. */}
                            <div className="panel-block" style={{paddingLeft: "20px"}}>
                                <div className="column is-full">
                                    <div className="field has-addons">
                                        <div className="control">
                                            <label className="radio">
                                                <input
                                                    type="radio"
                                                    name="inputType"
                                                    value="Fasta"
                                                    checked={selectedInputType === "Fasta"}
                                                    onChange={handleInputTypeChange}
                                                />
                                                <span style={{marginLeft: "5px"}}>
                                                    Fasta
                                                </span>
                                            </label>
                                        </div>
                                    </div>
                                    <div className="field has-addons">
                                        <div className="control">
                                            <label className="radio">
                                                <input
                                                    type="radio"
                                                    name="inputType"
                                                    value="Genbank"
                                                    checked={selectedInputType === "Genbank"}
                                                    onChange={handleInputTypeChange}
                                                />
                                                <span style={{marginLeft: "5px"}}>
                                                    Genbank
                                                </span>
                                            </label>
                                        </div>
                                    </div>
                                </div>
                            </div>

                            {/* Input block. */}
                            {selectedInputType === "Fasta" ? (
                                <FastaInputContainer src={src} setSrc={setSrc} />
                            ) : (
                                <GenbankInputContainer src={src} setSrc={setSrc} />
                            )}

                            {/* Submit button that sends input to API. */}
                            <div className="panel-block">
                                <div className="column is-full">
                                    <div className="field is-grouped is-grouped-right">
                                        <div className="control">
                                            {handleLoadExample && (
                                                <button
                                                    className="button is-secondary"
                                                    style={{marginRight: "10px"}}
                                                    onClick={() => {
                                                        handleLoadExample(setSrc);
                                                    }}
                                                >
                                                    Load Example
                                                </button>
                                            )}
                                            <button
                                                className="button is-secondary"
                                                style={{marginRight: "10px"}}
                                                onClick={() => {
                                                    handleRefresh();
                                                }}
                                            >
                                                Refresh
                                            </button>
                                            <button 
                                                className="button is-primary"
                                                onClick={() => {
                                                    if (src === "") {
                                                        toast.error("No input provided.");
                                                        return;
                                                    }
                                                    handleSubmit();
                                                }}
                                            >
                                                Submit
                                            </button>
                                        </div>
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

export default Input;