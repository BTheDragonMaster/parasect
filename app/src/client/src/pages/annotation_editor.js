import React, { useEffect, useState } from 'react';
import { useParams } from 'react-router-dom';
import { toast } from 'react-toastify';
import { Box, IconButton, Divider, Typography, Button, Modal, Tooltip } from '@mui/material';
import { MdClose } from 'react-icons/md';

import Loading from '../components/Loading';
import ProteinTile from '../components/ProteinTile';

import Turnstile from 'react-turnstile';


const SITE_KEY = process.env.REACT_APP_TURNSTILE_SITE_KEY;


function SubmitAnnotationsModal({ open, onClose, proteinAnnotations }) {
    const [captchaToken, setCaptchaToken] = useState(null);
    const [submitting, setSubmitting] = useState(false);
    const [turnstileKey, setTurnstileKey] = useState(0); // force reset when modal re-opens

    useEffect(() => {
        if (open) {
            // reset widget each time the modal opens
            setCaptchaToken(null);
            setTurnstileKey(k => k + 1);
        }
    }, [open]);

    {/* Submit updated protein data */}
    const handleSubmit = async () => {
        if (submitting) return; // prevent multiple submissions

        if (Object.keys(proteinAnnotations).length === 0) {
            toast.error("No annotations to submit");
            return;
        };

        if (!captchaToken) return;
        setSubmitting(true);

        try {
            const res = await fetch("/api/submit_annotations", {
                method: "POST",
                headers: {"Content-Type": "application/json"},
                body: JSON.stringify({
                    annotations: proteinAnnotations,
                    turnstileToken: captchaToken
                }),
                credentials: 'include' // include cookies
            });

            if (!res.ok) {
                const err = await res.json().catch(() => ({}));
                throw new Error(err?.error || "Submission failed")
            }

            // success UI
            onClose();
        } catch(e) {
            console.error(e);
            toast.error(`Error: ${e.message}`, { autoClose: false });
        } finally {
            setSubmitting(false);
            setCaptchaToken(null);
            setTurnstileKey(k => k + 1);  // reset widget after submit
        }
    };

    return (
        <Modal open={open} onClose={onClose}>
            <Box
                width={800}
                bgcolor='white.main'
                mx='auto'
                my={10}
                borderRadius={4}
                boxShadow={3}
            >
                <Box
                    sx={{
                        backgroundColor: 'secondary.main',
                        borderTopLeftRadius: '14px',
                        borderTopRightRadius: '14px',
                        display: 'flex',
                        justifyContent: 'space-between',
                        alignItems: 'center',
                        padding: 1,
                    }}
                >
                    <Typography 
                        variant='h5' 
                        gutterBottom
                        sx={{ 
                            color: 'black.main', 
                            textAlign: 'center',
                            pl: 2,
                            pt: 2, 
                        }}
                    >
                        Submit annotations
                    </Typography>
                    <IconButton onClick={onClose}>
                        <MdClose size={24} />
                    </IconButton>
                </Box>

                <Box
                    sx={{
                        display: 'flex',
                        flexDirection: 'column',
                        alignItems: 'left',
                        gap: 3,
                        p: 4,
                    }}
                >
                    {/* Turnstile widget */}
                    <Turnstile
                        key={turnstileKey}
                        sitekey={SITE_KEY}
                        onVerify={(token) => setCaptchaToken(token)}
                        onExpire={() => setCaptchaToken(null)}
                        onError={() => setCaptchaToken(null)}
                        // options={{ theme: 'auto', appearance: 'always' }} // optional
                    />
                    <Button
                        variant="contained"
                        color="primary"
                        onClick={handleSubmit}
                        disabled={!captchaToken || submitting}
                    >
                        {submitting ? 'Submitting…' : 'Submit'}
                    </Button>
                </Box>
            </Box>
        </Modal>
    )  
};


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

    const [proteinAnnotations, setProteinAnnotations] = useState({});

    // submission modal state
    const [openAnnotationsSubmissionModal, setOpenAnnotationsSubmissionModal] = useState(false);

    const handleOpenAnnotationsSubmissionModal = () => {
        setOpenAnnotationsSubmissionModal(true);
    };

    const handleCloseAnnotationsSubmissionModal = () => {
        setOpenAnnotationsSubmissionModal(false);
    };

    {/* For collecting protein annotations */}

    const handleProteinAnnotationChange = (proteinId, data) => {
    setProteinAnnotations((prev) => {
        const updated = { ...prev };

        // Check if domains object exists and has keys
        const hasDomainAnnotations = data.domains && Object.keys(data.domains).length > 0;

        if (data && hasDomainAnnotations) {
            updated[proteinId] = data;
        } else {
            // Remove if no domain annotations (ignore synonym)
            delete updated[proteinId];
        }

        return updated;
    });
};

    const OpenAnnotationsSubmissionsModal = () => {
        setOpenAnnotationsSubmissionModal(true);
    };

    // fetch results from local storage
    useEffect(() => {
        let intervalId;

        const fetchResult = async () => {
            try {
                const response = await fetch(`/api/retrieve/${jobId}`);
                if (!response.ok) {
                    throw new Error('failed to fetch results');
                }

                const data = await response.json();

                if (data.status === 'success') {
                    const results = data['payload']['results']

                    // set states
                    setResults(results);
                    setIsLoading(false);
                    clearInterval(intervalId);
                } else if (data.status === 'failure') {
                    throw new Error(data.message);
                } // else keep polling
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
            }
        };

        if (jobId) {
            // poll every second (1000 milliseconds)
            intervalId = setInterval(fetchResult, 1000);
        }

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
                <p>Extracting A-domains and making PARAS predictions...</p>
            </Box>
        );
    }

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
    }

    // render results if available
    return (
        <>
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
                    <Box sx={{ display: 'flex', flexDirection: 'row', gap: 1 }}>
                        <Typography variant='h4' gutterBottom>
                            Extracted adenylation domains
                        </Typography>
                    </Box>
                    <Divider />

                    <Box sx={{ mt: 4 }}>
                        <Typography variant='body1' gutterBottom>
                            In total, {results.length} of the submitted proteins contain adenylation domains. The tiles below show the domains for each protein. You can scroll horizontally to view all proteins.
                        </Typography>

                        <Typography variant='body1' gutterBottom>
                            Domains which already exist in the PARAS/PARASECT dataset are displayed in grey. New domains are displayed in yellow. Please review the new domains and provide annotations where possible. You can proceed with submitting domains once they are annotated.
                        </Typography>

                        <Typography variant='body1' gutterBottom>
                        </Typography>
                        <Box sx={{ mt: 4 }}>
                            <Tooltip title={Object.keys(proteinAnnotations).length === 0 ? "Nothing to submit! Have you annotated all new domains?" : ""} arrow>
                                {/* Wrap Button in a span to support tooltip for disabled state */}
                                <span>
                                    <Button
                                        variant="contained"
                                        color="primary"
                                        onClick={OpenAnnotationsSubmissionsModal}
                                        disabled={Object.keys(proteinAnnotations).length === 0}
                                    >
                                        Proceed with submission
                                    </Button>
                                </span>
                            </Tooltip>
                        </Box>
                    </Box>
                </Box>

                {/* display results in a row, one item per protein */}
                <Box
                    sx={{
                        overflowY: 'auto',
                        overflowX: 'auto',
                        backgroundColor: 'white.main',
                        flexDirection: 'row',
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
                        <ProteinTile
                            key={index}
                            proteinResult={result}
                            onUpdateAnnotation={handleProteinAnnotationChange} />
                    ))}
                </Box>
            </Box>

            {/* Dialogue modal */}
            <SubmitAnnotationsModal 
                open={openAnnotationsSubmissionModal}
                onClose={handleCloseAnnotationsSubmissionModal}
                proteinAnnotations={proteinAnnotations}
            />
        </>
        
    );
};

export default AnnotationEditor;