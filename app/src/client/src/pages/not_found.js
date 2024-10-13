import React from 'react';
import { Link } from 'react-router-dom';
import { Box, Typography, Button } from '@mui/material';

const NotFound = () => {
    return (
        <Box
            display='flex'
            flexDirection='column'
            justifyContent='center'
            alignItems='center'
            minHeight='80vh' // vertically centers content without taking full screen
            padding={4} // arger padding for better spacing
        >
            <Box
                display='flex'
                flexDirection='column'
                justifyContent='center'
                alignItems='center'
                mb={2}
            >
                <Typography variant='h4' component='div' fontWeight='bold' gutterBottom>
                    404 Not Found
                </Typography>
                <Box
                    component='img'
                    src='/paras_error.png'
                    alt='PARAS Error'
                    sx={{ width: 300 }}
                />
            </Box>
            <Button component={Link} to='/' variant='contained' color='primary'>
                Go back to the home page
            </Button>
        </Box>
    );
};

export default NotFound;