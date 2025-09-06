import React from 'react';
import { Box, Typography, Button, Link, Tooltip } from '@mui/material';
import ExitIcon from '@mui/icons-material/ExitToApp';

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
                    sx={{ maxWidth: 300 }}
                />
                <Box
                    component='img'
                    src='/parasect.png'
                    alt='PARASECT logo'
                    sx={{ maxWidth: 300 }}
                />
            </Box>

            <Typography variant='h3' component='div' fontWeight='bold' gutterBottom>
                Welcome to PARAS(ECT)!
            </Typography>

            <Typography variant='h6' component='div' color='textSecondary' align='center' gutterBottom>
                Discover our adenylation domain substrate specificity prediction models.
            </Typography>

            <Box
                sx={{
                    mt: 4,
                    mb: 2,
                    textAlign: 'center',
                    display: 'flex',
                    flexDirection: 'row',
                    gap: 2
                }}
            >
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
                <Button
                    variant='contained'
                    color='primary'
                    size='large'
                    href='/retrieve'
                    sx={{ marginBottom: 3 }}
                >
                    <Typography sx={{ color: 'white.main' }}>
                        Retrieve results
                    </Typography>
                </Button>
                <Button
                    variant='contained'
                    color='primary'
                    size='large'
                    href='/data_annotation'
                    sx={{ marginBottom: 3 }}
                >
                    <Typography sx={{ color: 'white.main' }}>
                        Annotate domain
                    </Typography>
                </Button>
            </Box>

            <Typography variant='body1' align='center' color='textSecondary' gutterBottom>
                <Box sx={{ display: 'flex', alignItems: 'center' }}>
                    Want to learn more about the research behind PARAS and PARASECT?
                    <Tooltip title='Opens new tab to bioRxiv.' arrow>
                        <Link 
                            href='https://www.biorxiv.org/content/10.1101/2025.01.08.631717v1'
                            underline='hover' 
                            target='_blank'
                            sx={{ marginLeft: 1, fontWeight: 'bold' }}
                        >
                            <Box sx={{ display: 'flex', alignItems: 'center', gap: 0.5 }}>
                                Read our publication 
                                <ExitIcon fontSize='small' />
                            </Box>
                        </Link>
                    </Tooltip>
                </Box>
            </Typography>

            <Typography variant='body1' align='center' color='textSecondary' gutterBottom>
                <Box sx={{ display: 'flex', alignItems: 'center' }}>
                    Did you find PARAS or PARASECT useful?
                    <Tooltip title='Opens new tab to bioRxiv.' arrow>
                        <Link 
                            href='https://www.biorxiv.org/content/10.1101/2025.01.08.631717v1'
                            underline='hover' 
                            target='_blank'
                            sx={{ marginLeft: 1, fontWeight: 'bold' }}
                        >
                            <Box sx={{ display: 'flex', alignItems: 'center', gap: 0.5 }}>
                                Please cite our publication
                                <ExitIcon fontSize='small' />
                            </Box>
                        </Link>
                    </Tooltip>
                </Box>
            </Typography>

            <Typography variant='body1' align='center' color='textSecondary'>
                <Box sx={{ display: 'flex', alignItems: 'center' }}>
                    Have something to contribute?
                    <Tooltip title='Opens new tab to GitHub.' arrow>
                        <Link 
                            href='https://github.com/bthedragonmaster/parasect'
                            target='_blank'
                            underline='hover'
                            sx={{ marginLeft: 1, fontWeight: 'bold' }}
                        >
                            <Box sx={{ display: 'flex', alignItems: 'center', gap: 0.5 }}>
                                Visit our GitHub page
                                <ExitIcon fontSize='small' />
                            </Box>
                        </Link>
                    </Tooltip>
                </Box>
            </Typography>
        </Box>
    );
};

export default Home;
