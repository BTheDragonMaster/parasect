import React, { useEffect, useMemo, useState } from 'react';
import { useParams } from 'react-router-dom';
import { toast } from 'react-toastify';
import { Box, IconButton, Divider, Typography, FormControl, InputLabel, Select, MenuItem } from '@mui/material';
import { FaDownload, FaCopy } from 'react-icons/fa';

import Loading from '../components/Loading';
import ResultTile from '../components/ResultTile';

/**
 * Component to display the results of the prediction.
 * 
 * @returns {React.ReactElement} - The results component.
 */
const Results = () => {
    // get job ID from URL
    const { jobId } = useParams();

    // state to store results
    const [results, setResults] = useState(null);

    // state to keep track of loading state
    const [isLoading, setIsLoading] = useState(true);

    // what the server last said about the job: null until the first response
    // lands, then 'pending' while it is genuinely still running. This page
    // serves both flows (a fresh submit and a retrieve of an older job) and
    // nothing in the route tells them apart, so the status is what decides
    // which loading message is honest.
    const [jobStatus, setJobStatus] = useState(null);

    // how to order the protein/gene groups of prediction cards
    // ('default' = order as returned by the server and for submit_quick/submit_domain
    // this reflects submission order; for a FASTA/GBK upload it's the order in
    // which the server extracted the domains, not necessarily file/genomic order)
    const [sortBy, setSortBy] = useState('default');

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
                setJobStatus(data.status);

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
            // fetch straight away, then poll every second (1000 milliseconds).
            // Without the leading call an already-finished job still sits on
            // the loading screen for a full second before its results appear.
            fetchResult();
            intervalId = setInterval(fetchResult, 1000);
        };

        // clear interval when component unmounts
        return () => clearInterval(intervalId);

    }, [jobId]);

    // group prediction cards by protein/gene, then order those groups according
    // to sortBy. Domains within a gene are always ordered by domain_nr (their
    // position along that protein). Any other order within a gene would be
    // confusing, so that part isn't user-configurable.
    const groups = useMemo(() => {
        if (!results) return [];

        const groupOrder = []; // order of first appearance, for the 'default' sort
        const byProtein = new Map();
        results.forEach((result) => {
            const name = result['domain_name'];
            if (!byProtein.has(name)) {
                byProtein.set(name, []);
                groupOrder.push(name);
            }
            byProtein.get(name).push(result);
        });

        const topConfidence = (name) => Math.max(
            ...byProtein.get(name).map((r) => parseFloat(r['predictions']?.[0]?.probability ?? 0))
        );

        // earliest position along the DNA among a gene's domains, from a GBK
        // upload (see domain_genomic_position); Infinity when unknown, so
        // genes without a known position sort to the end rather than the start
        const earliestGenomicPosition = (name) => Math.min(
            ...byProtein.get(name).map((r) => (
                r['domain_genomic_position'] === null || r['domain_genomic_position'] === undefined
                    ? Infinity
                    : r['domain_genomic_position']
            ))
        );

        let proteinNames;
        switch (sortBy) {
            case 'protein_asc':
                proteinNames = [...groupOrder].sort((a, b) => a.localeCompare(b));
                break;
            case 'protein_desc':
                proteinNames = [...groupOrder].sort((a, b) => b.localeCompare(a));
                break;
            case 'confidence_desc':
                proteinNames = [...groupOrder].sort((a, b) => topConfidence(b) - topConfidence(a));
                break;
            case 'domain_count_desc':
                proteinNames = [...groupOrder].sort((a, b) => byProtein.get(b).length - byProtein.get(a).length);
                break;
            case 'genomic_asc':
                proteinNames = [...groupOrder].sort((a, b) => earliestGenomicPosition(a) - earliestGenomicPosition(b));
                break;
            default:
                proteinNames = groupOrder;
        }

        return proteinNames.map((name) => ({
            proteinName: name,
            items: [...byProtein.get(name)].sort((a, b) => a['domain_nr'] - b['domain_nr']),
        }));
    }, [results, sortBy]);

    // the "order in DNA" option only makes sense (and is only available) for
    // GBK uploads. The server doesn't track genomic position for plain
    // FASTA/protein input or signatures submitted directly
    const hasGenomicPositions = useMemo(
        () => !!results?.some((r) => r['domain_genomic_position'] !== null && r['domain_genomic_position'] !== undefined),
        [results]
    );

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
                <p>{jobStatus === 'pending' ? 'Making predictions...' : 'Retrieving predictions...'}</p>
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
                <Box sx={{ mt: 1, mb: 1 }}>
                    <FormControl size='small' sx={{ minWidth: 240 }}>
                        <InputLabel id='sort-by-label'>Group order</InputLabel>
                        <Select
                            labelId='sort-by-label'
                            label='Group order'
                            value={sortBy}
                            onChange={(e) => setSortBy(e.target.value)}
                        >
                            <MenuItem value='default'>Default (server order)</MenuItem>
                            {hasGenomicPositions && (
                                <MenuItem value='genomic_asc'>Order in DNA (GBK only)</MenuItem>
                            )}
                            <MenuItem value='protein_asc'>Protein name (A-Z)</MenuItem>
                            <MenuItem value='protein_desc'>Protein name (Z-A)</MenuItem>
                            <MenuItem value='confidence_desc'>Highest confidence first</MenuItem>
                            <MenuItem value='domain_count_desc'>Most domains first</MenuItem>
                        </Select>
                    </FormControl>
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
                        In total, {results.length} prediction(s) were made across {groups.length} protein(s). Domains from the same protein are grouped together and shaded the same color, in their order along that protein. You can scroll horizontally to view all predictions, and use the dropdown above to change the order proteins are grouped in.
                    </Typography>
                    <Typography variant='body1' gutterBottom>
                        You can download the results as a JSON file using the download button above next to the header.
                    </Typography>
                    <Typography variant='body1' gutterBottom>
                        You can use the job ID to retrieve the results at a later time. All jobs are automatically deleted after 7 days.
                    </Typography>
                </Box>
            </Box>



            {/* display results in a row, one group of cards per protein/gene */}
            <Box
                sx={{
                    overflowY: 'auto',
                    overflowX: 'auto',
                    backgroundColor: 'background.default',
                    display: 'flex',
                    alignItems: 'flex-start',
                    gap: '32px',
                    paddingLeft: '30px',
                    paddingRight: '30px',
                    paddingBottom: '20px',

                    // always show scrollbar
                    '&::-webkit-scrollbar': {
                        display: 'block',
                    },

                    // scrollbar style
                    '&::-webkit-scrollbar-thumb': {
                        backgroundColor: 'surface.borderStrong',
                        borderRadius: '10px',
                    },
                }}
            >
                {groups.map((group, groupIndex) => (
                    <Box
                        key={group.proteinName}
                        sx={{ display: 'flex', flexDirection: 'column', flexShrink: 0 }}
                    >
                        <Box
                            sx={{
                                display: 'flex',
                                alignItems: 'baseline',
                                gap: 1,
                                mb: 1,
                                px: 1,
                                py: 0.5,
                                borderRadius: '8px 8px 0 0',
                                // the gene's colour lives in a solid bar, not in a
                                // tint of the panel: a 7%-opacity fill was all but
                                // invisible against the page
                                borderBottom: '3px solid',
                                borderColor: (theme) => theme.palette.geneBands[groupIndex % 2].rail,
                            }}
                        >
                            <Typography variant='subtitle2' sx={{ fontWeight: 700 }}>
                                {group.proteinName}
                            </Typography>
                            <Typography variant='caption' color='textSecondary'>
                                {group.items.length} domain{group.items.length > 1 ? 's' : ''}
                            </Typography>
                        </Box>
                        <Box
                            sx={{
                                display: 'flex',
                                gap: '20px',
                                p: 1.5,
                                borderRadius: '0 14px 14px 14px',
                                backgroundColor: (theme) => theme.palette.geneBands[groupIndex % 2].surface,
                                borderLeft: '4px solid',
                                borderColor: (theme) => theme.palette.geneBands[groupIndex % 2].rail,
                            }}
                        >
                            {group.items.map((result) => (
                                <ResultTile key={`${group.proteinName}-${result['domain_nr']}`} result={result} />
                            ))}
                        </Box>
                    </Box>
                ))}
            </Box>
        </Box>
    );
};

export default Results;