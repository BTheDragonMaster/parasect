import React, { useState } from "react";

const FastaInputContainer = () => {
    const [src, setSrc] = useState("");

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
                            placeholder="Enter Fasta sequence here or upload a file..."
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

const GenbankInputContainer = () => {
    const [src, setSrc] = useState("");

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

const Input = () => {
    const [selectedInputType, setSelectedInputType] = useState("Fasta");

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
                    </div>
                </div>

                {/* Add a radio button to select the input type: Fasta or Genbank. */}
                <div className="panel-block" style={{paddingLeft: "20px"}}>
                    <div className="column is-full">
                        <div className="field has-addons">
                            <div className="control">
                                <label class="radio">
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
                                <label class="radio">
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
                    <FastaInputContainer />
                ) : (
                    <GenbankInputContainer />
                )}
            </div>
        </div>
    );
};

export default Input;