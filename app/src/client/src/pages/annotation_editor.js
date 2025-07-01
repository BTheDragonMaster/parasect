import React, { useEffect, useState } from 'react';
import { useParams } from 'react-router-dom';
import { toast } from 'react-toastify';
import { Box, IconButton, Divider, Typography } from '@mui/material';
import { FaDownload, FaCopy } from 'react-icons/fa';

import Loading from '../components/Loading';
import ResultTile from '../components/ResultTile';

/**
 * Component to display the results of the prediction.
 *
 * @returns {React.ReactElement} - The results component.
 */
const AnnotationEditor = () => {
    // get job ID from URL
    const { jobId } = useParams();

    // state to store results
    const [results, setResults] = useState(null);

    // state to keep track of loading state
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
                    const results = data['payload']['results']

                    // sort all predictions by probability
                    results.forEach(result => {
                        result['predictions'].sort((a, b) => b['probability'] - a['probability']);
                    });

                    // round all probabilities to 2 decimal places
                    results.forEach(result => {
                        result['predictions'].forEach(prediction => {
                            prediction['probability'] = prediction['probability'].toFixed(3);
                        });
                    });

                    // set states
                    setResults(results);
                    setIsLoading(false);
                    clearInterval(intervalId);
                } else if (data.status === 'failure') {
                    throw new Error(data.message);
                }; // else keep polling
            } catch (error) {
                toast.error(
                    <>
                      {error.message}<br /><br />
                      If you feel this is an error, or if you need assistance, please contact the developers in GitHub issues by selecting 'Report an issue' in the app bar at the top left of this page and posting your issue or question.
                    </>,
                    { autoClose: false }
                );
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

    // render loading spinner while fetching results
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

    // render message if no results are found
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

    // render results if available
    return (
        <Box
            display='flex'
            flexDirection='column'
            overflow='hidden'
        >
            <Box
                display='flex'
                flexDirection='column'
                alignItems='left'
                margin={4}
            >

                {/* header with job ID and download button */}
                <Box sx={{ display: 'flex', flexDirection: 'row', gap: 1 }}>
                    <Typography variant='h4' gutterBottom>
                        Results
                    </Typography>
                    <Box>

                        {/* download button */}
                        <IconButton
                            onClick={() => {
                                const blob = new Blob([JSON.stringify(results)], { type: 'application/json' });
                                const url = URL.createObjectURL(blob);
                                const a = document.createElement('a');
                                a.href = url;
                                a.download = 'results.json';
                                a.click();
                            }
                        }>
                            <FaDownload />
                        </IconButton>

                    </Box>
                </Box>
                <Typography variant='body1' gutterBottom>
                    <IconButton
                        onClick={() => {
                            navigator.clipboard.writeText(jobId);
                            toast.success('Copied the job ID to clipboard!');
                        }}
                    >
                        <FaCopy size={15} style={{ paddingBottom: '3px' }} />
                    </IconButton>
                    {`Job ID: ${jobId}`}
                </Typography>
                <Divider />

                <Box sx={{ mt: 4 }}>
                    <Typography variant='body1' gutterBottom>
                        In total, {results.length} prediction(s) were made. The tiles below show the predictions for each domain. You can scroll horizontally to view all predictions.
                    </Typography>
                    <Typography variant='body1' gutterBottom>
                        You can download the results as a JSON file using the download button above next to the header.
                    </Typography>
                    <Typography variant='body1' gutterBottom>
                        You can use the job ID to retrieve the results at a later time. All jobs are automatically deleted after 7 days.
                    </Typography>
                </Box>
            </Box>



            {/* display results in a row, one item per domain with prediction */}
            <Box
                sx={{
                    overflowY: 'auto',
                    overflowX: 'auto',
                    backgroundColor: 'white.main',
                    flexDirection: 'column',
                    display: 'flex',
                    gap: '20px',
                    paddingLeft: '30px',
                    paddingRight: '30px',
                    paddingBottom: '20px',

                    // always show scrollbar
                    '&::-webkit-scrollbar': {
                        display: 'block',
                    },

                    // scrollbar style
                    '&::-webkit-scrollbar-thumb': {
                        backgroundColor: '#ceccca',
                        borderRadius: '10px',
                    },
                }}
            >
                {results.map((result, index) => (
                    <ResultTile key={index} result={result} />
                ))}
            </Box>
        </Box>
    );
};

export default AnnotationEditor;