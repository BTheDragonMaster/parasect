import React from 'react';
import { Box, Typography, Button, Link } from '@mui/material';

/**
 * Home component that displays the home page content.
 * 
 * @returns {React.ReactElement} - The component showing the home page content.
 */
const Home = () => {
    return (
        <Box
            display='flex'
            flexDirection='column'
            justifyContent='center'
            alignItems='center'
            minHeight='80vh' // vertically centers content without taking full screen
            padding={4} // larger padding for better spacing
        >
            <Box display='flex' flexDirection='row' justifyContent='center' alignItems='center'>
                <Box
                    component='img'
                    src='/paras.png'
                    alt='PARAS logo'
                    sx={{ width: 300 }}
                />
                <Box
                    component='img'
                    src='/parasect.png'
                    alt='PARASECT logo'
                    sx={{ width: 300 }}
                />
            </Box>

            <Typography variant='h3' component='div' fontWeight='bold' gutterBottom>
                Welcome to PARAS(ECT)!
            </Typography>

            <Typography variant='h6' component='div' color='textSecondary' align='center' gutterBottom>
                Discover our adenylation domain substrate specificity prediction models.
            </Typography>

            <Box mt={4} mb={2} textAlign='center'>
                <Button
                    variant='contained'
                    color='primary'
                    size='large'
                    href='/submit'
                    sx={{ marginBottom: 3 }}
                >
                    <Typography sx={{ color: 'white.main' }}>
                        Start predicting
                    </Typography>
                </Button>
            </Box>

            <Typography variant='body1' align='center' color='textSecondary' gutterBottom>
                Want to learn more about the research behind PARAS and PARASECT?
                <Link 
                    href='/publication' 
                    underline='hover' 
                    sx={{ marginLeft: 1, fontWeight: 'bold' }}
                >
                    Read our publication.
                </Link>
            </Typography>

            <Typography variant='body1' align='center' color='textSecondary'>
                Have something to contribute?
                <Link 
                    href='https://github.com/bthedragonmaster/parasect'
                    target='_blank'
                    underline='hover'
                    sx={{ marginLeft: 1, fontWeight: 'bold' }}
                >
                    Visit our GitHub page.
                </Link>
            </Typography>
        </Box>
    );
};

export default Home;
