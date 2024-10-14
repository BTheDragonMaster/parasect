import React, { useEffect, useState } from 'react';
import { useParams } from 'react-router-dom';
import { toast } from 'react-toastify';
import { Box } from '@mui/material';
import Loading from '../components/Loading';

const Results = () => {
    const { jobId } = useParams();
    const [results, setResults] = useState(null);
    const [isLoading, setIsLoading] = useState(true);

    // fetch results from local storage
    useEffect(() => {
        let intervalId;

        const fetchResult = async () => {
            try {
                const response = await fetch(`/api/retrieve/${jobId}`);
                if (!response.ok) {
                    throw new Error('failed to fetch results');
                };

                const data = await response.json();

                if (data.status === 'success') {
                    setResults(data['payload']['results']);
                    setIsLoading(false);
                    clearInterval(intervalId);
                } else if (data.status === 'failure') {
                    throw new Error(data.message);
                }; // else keep polling
            } catch (error) {
                toast.error(error.message);
                setIsLoading(false);
                clearInterval(intervalId);
            };
        };

        if (jobId) {
            // poll every second (1000 milliseconds)
            intervalId = setInterval(fetchResult, 1000);
        };

        // clear interval when component unmounts
        return () => clearInterval(intervalId);

    }, [jobId]);

    if (isLoading) {
        return (
            <Box
                display='flex'
                flexDirection='column'
                justifyContent='center'
                alignItems='center'
                minHeight='80vh'
            >
                <Loading 
                    frame1='paras_loading_1.png' 
                    frame2='paras_loading_2.png' 
                />
                <p>Making predictions...</p>
            </Box>
        );  
    };

    if (!results) {
        return (
            <Box
                display='flex'
                justifyContent='center'
                alignItems='center'
                minHeight='80vh'
            >
                <p>No results found for job ID {jobId}</p>
            </Box>
        );
    };

    return (
        <Box>
            <pre>{JSON.stringify(results, null, 2)}</pre>
        </Box>
    );
};

export default Results;