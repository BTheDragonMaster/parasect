import React, { useEffect } from "react";
import { toast } from "react-toastify";
import { Box } from "@mui/material";

const Results = ({ results, setResults }) => {

    // fetch results from local storage
    useEffect(() => {
        const storedResults = localStorage.getItem("results");
        if (storedResults) {
            setResults(JSON.parse(storedResults));
        } else {
            toast.error("No results found!");
        }
    }, []);

    return (
        <Box>
            <pre>{JSON.stringify(results, null, 2)}</pre>
        </Box>
    );
};

export default Results;