import React, { useState } from "react";
import DropdownList from "react-widgets/DropdownList";

const Results = ({ results }) => {
    const domain_ids = results.map((result) => result.domain_id);
    const [selectedDomainId, setSelectedDomainId] = useState(domain_ids[0]);

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
                    <button
                        className="button is-primary"
                        style={{ marginBottom: "10px" }}
                        onClick={handleDownload}
                    >
                        Download results
                    </button>
                    <DropdownList
                        data={domain_ids}
                        defaultValue={selectedDomainId}
                        onChange={(value) => setSelectedDomainId(value)}
                    />
                    <div>
                        {results.map((result, result_index) => {
                            if (result.domain_id === selectedDomainId) {
                                return (
                                    <div 
                                        key={result_index}
                                        style={{marginTop: "10px"}}
                                    >
                                        <div style={{ overflowX: "auto", overflowY: "auto" }}>
                                            <table style={{ borderCollapse: "collapse" }}>
                                                <thead>
                                                    <tr>
                                                        <th />
                                                        <th style={{textAlign: "left", whiteSpace: "nowrap"}}>Name</th>
                                                        <th style={{textAlign: "left", whiteSpace: "nowrap"}}>Score</th>
                                                    </tr>
                                                </thead>
                                                <tbody>
                                                    {result.data.predictions.map((item, item_index) => (
                                                        <tr key={item_index}>
                                                            <td>{item_index + 1}</td>
                                                            <td style={{textAlign: "left", whiteSpace: "nowrap"}}>{item[1]}</td>
                                                            <td style={{textAlign: "left", whiteSpace: "nowrap"}}>{item[0].toFixed(3)}</td>
                                                        </tr>
                                                    ))}
                                                </tbody>
                                            </table>
                                        </div>
                                    </div>
                                );
                            }
                        })}
                    </div>
                </div>
            )}
        </div>
    );
};

export default Results;