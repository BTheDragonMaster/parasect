import React, { useState, useEffect } from "react";
import DropdownList from "react-widgets/DropdownList";
import Plot from "react-plotly.js";

const layout = {
    autoresize: true,
    autosize: true,
    xaxis: {
        title: {
            text: "Prediction value",
            standoff: 20,
            font: { size: 16, family: "Arial" }
        },
        autorange: "reversed",
        range: [0, 1.1],
        tickformat: ",.1",
        hoverformat: ",.3",
        automargin: true,
        titlefont: { size: 16, family: "Arial" },
        tickfont: { size: 16, family: "Arial" },
        ticklen: 10,
        side: "top",
    },
    yaxis: {
        automargin: true,
        titlefont: { size: 16, family: "Arial" },
        tickfont: { size: 16, family: "Arial" },
        ticklen: 10,
        side: "right"
    },
    margin: { l: 0, r: 0, b: 50, t: 0, pad: 0 },
};

function json2tsv(results) {
    let tsv = "domain_id\t";
    let predictions = results[0].data.predictions;
    for (let i = 0; i < predictions.length; i++) {
        tsv += `prediction_${i+1}\tvalue_${i+1}\t`;
    }
    tsv += "\n";
    for (let result of results) {
        tsv += `${result.domain_id}\t`;
        for (let prediction of result.data.predictions) {
            tsv += `${prediction[1]}\t${prediction[0]}\t`;
        }
        tsv += "\n";
    }
    return tsv;
};

const Results = ({ results }) => {
    const domain_ids = results.map((result) => result.domain_id);
    const [selectedDomainId, setSelectedDomainId] = useState(domain_ids[0]);

    const [probabilities, setProbabilities] = useState([]);
    const [labels, setLabels] = useState([]);

    useEffect(() => {
        if (results) {
            const selectedResult = results.find((result) => result.domain_id === selectedDomainId);
            const predictions = selectedResult.data.predictions;
            const probabilities = predictions.map((item) => item[0]);
            const labels = predictions.map((item) => item[1]);
            probabilities.reverse();
            labels.reverse();
            setProbabilities(probabilities);
            setLabels(labels);
        }
    }, [results, selectedDomainId]);

    const handleDownload = () => {
        const filename = "results.json";
        const data = JSON.stringify(results, null, 4);
        const blob = new Blob([data], { type: "application/json" });
        const url = URL.createObjectURL(blob);
        const link = document.createElement("a");
        link.href = url;
        link.download = filename;
        link.click();
    };

    

    return (
        <div>
            {results && (
                <div>
                    <div style={{ margin: "15px" }}>
                        {`
                            There are ${results.length} domain(s) found in the 
                            input for which substrate specificity predictions are 
                            available. The predictions are based on the trained 
                            model and the results are shown below. You can select 
                            a domain from the dropdown list to view the predictions.                            
                        `}
                    </div>
                    <div style={{ margin: "15px" }}>
                        {`
                            If you have selected any additional options for metadata
                            in the Options window before submitting the input, you can
                            download this metadata along with the predictions as JSON.
                            The predictions alone can also be downloaded as a tab-separated
                            values (TSV) file.
                        `}
                    </div>
                <div 
                    style={{
                        display: "flex",
                        flexDirection: "row",
                        justifyContent: "center",
                        alignItems: "left",
                        width: "100%",
                        height: "100%",
                        padding: "10px",
                        gap: "30px",
                    }}
                >
                    <div
                        style={{
                            display: "flex",
                            flexDirection: "column",
                            alignItems: "left",
                            width: "70%",
                            height: "100%",
                        }}
                    >
                        <DropdownList
                            data={domain_ids}
                            defaultValue={selectedDomainId}
                            onChange={(value) => setSelectedDomainId(value)}
                            style={{ marginBottom: "10px" }}
                        />
                        <Plot
                            data={[
                                {
                                    x: probabilities,
                                    y: labels,
                                    type: "bar",
                                    orientation: "h",
                                    marker: { color: "#f28732" }
                                }
                            ]}
                            layout={layout}
                            config={{ responsive: true }}
                        />
                    </div>
                    <div 
                        style={{
                            display: "flex",
                            flexDirection: "column",
                            justifyContent: "left",
                            alignItems: "left",
                            width: "30%",
                            height: "100%",
                        }}
                    >
                        <button
                            className="button"
                            style={{ 
                                backgroundColor: "#ccc",
                                color: "#000",
                                marginBottom: "10px",
                            }}
                            onClick={handleDownload}
                        >
                            Download all results (JSON)
                        </button>
                        <button
                            className="button"
                            style={{ 
                                backgroundColor: "#ccc",
                                color: "#000",
                                marginBottom: "10px",
                            }}
                            onClick={() => {
                                const tsv = json2tsv(results);
                                const blob = new Blob([tsv], { type: "text/tab-separated-values" });
                                const url = URL.createObjectURL(blob);
                                const link = document.createElement("a");
                                link.href = url;
                                link.download = "results.tsv";
                                link.click();
                            }}
                        >
                            Download predictions (TSV)
                        </button>
                    </div>
                </div>
                </div>
            )}
        </div>
    );
};

export default Results;