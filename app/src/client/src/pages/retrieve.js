import React, { useState } from 'react';
import { toast } from 'react-toastify';
import { Box, Typography, Divider, TextField, Button } from '@mui/material';

/**
 * Retrieve component that displays the retrieve page content.
 * 
 * @returns {React.ReactElement} - The component showing the retrieve page content.
 */
const Retrieve = () => {
    const [jobId, setJobId] = useState('');

    const handleRetrieve = () => {
        if (jobId) {
            window.location.href = `/results/${jobId}`;
        } else {
            toast.error('Please enter a Job ID!');
        }
    };

    return (
        <Box 
            display='flex' 
            flexDirection='column' 
            alignItems='left' 
            padding={4} 
            margin='auto'
        >
            
            <Typography variant='h4' gutterBottom>
                Retrieve results
            </Typography>
            <Divider />

            <Box sx={{ mt: 4 }}>
                <Typography variant='body1' gutterBottom>
                    All jobs are automatically deleted after 7 days. Please enter your Job ID below to retrieve the results.
                </Typography>
            </Box>

            {/* job ID input */}
            <TextField
                label="Job ID"
                variant="outlined"
                fullWidth
                value={jobId}
                onChange={(e) => setJobId(e.target.value)}
                margin="normal"
            />

            <Box mt={1} display='flex' justifyContent='left' width='100%' gap={2}>
                <Button 
                    variant="contained" 
                    color="primary" 
                    onClick={handleRetrieve}
                >
                    Submit
                </Button>
            </Box>
        </Box>
    );
};

export default Retrieve;
