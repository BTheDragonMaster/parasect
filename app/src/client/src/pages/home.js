import React from 'react';
import { Link as RouterLink } from 'react-router-dom';
import { Box, Typography, Link, Tooltip, Card, CardActionArea, Divider } from '@mui/material';
import ExitIcon from '@mui/icons-material/ExitToApp';
import UploadIcon from '@mui/icons-material/Upload';
import RetrieveIcon from '@mui/icons-material/GetApp';
import DatasetIcon from '@mui/icons-material/Dataset';
import QueryStatsIcon from '@mui/icons-material/QueryStats';
import HubIcon from '@mui/icons-material/Hub';
import ArrowForwardIcon from '@mui/icons-material/ArrowForward';

const links = [
    {
        text: 'Want to learn more about the research behind PARAS and PARASECT?',
        cta: 'Read our publication',
        href: 'https://pubs.acs.org/jaaucr/article/6/4/2315/5145284/PARAS-High-Accuracy-Machine-Learning-of-Substrate',
        tooltip: 'Opens new tab to JACS Au.',
    },
    {
        text: 'Did you find PARAS or PARASECT useful?',
        cta: 'Please cite our publication',
        href: 'https://pubs.acs.org/jaaucr/article/6/4/2315/5145284/PARAS-High-Accuracy-Machine-Learning-of-Substrate',
        tooltip: 'Opens new tab to JACS Au.',
    },
    {
        text: 'Have something to contribute?',
        cta: 'Visit our GitHub page',
        href: 'https://github.com/bthedragonmaster/parasect',
        tooltip: 'Opens new tab to GitHub.',
    },
];

const primaryAction = {
    to: '/submit',
    Icon: UploadIcon,
    label: 'Start predicting',
    blurb: 'Upload protein sequences or a GenBank file and predict what each adenylation domain activates.',
};

const secondaryActions = [
    {
        to: '/retrieve',
        Icon: RetrieveIcon,
        label: 'Retrieve results',
        blurb: 'Reopen an earlier run with its job ID. Jobs are kept for seven days.',
    },
    {
        to: '/query_database',
        Icon: QueryStatsIcon,
        label: 'Query database',
        blurb: 'Search the reference dataset by substrate, domain or taxonomy, and download what you find.',
    },
    {
        to: '/data_annotation',
        Icon: DatasetIcon,
        label: 'Annotate domains',
        blurb: 'Extract domains from your own sequences and contribute substrate annotations back.',
    },
    {
        to: '/network',
        Icon: HubIcon,
        label: 'Explore the network',
        blurb: 'Browse the reference domains as a similarity network, clustered by their signatures.',
    },
];

const LOGO_PLATE = {
    width: { xs: 116, sm: 150 },
    height: 'auto',
    borderRadius: 3,
    // give logos a white plate rather than letting it float on a dark page
    backgroundColor: '#FFFFFF',
    p: 1,
    border: '1px solid',
    borderColor: 'divider',
};

/**
 * Home component that displays the home page content.
 *
 * @returns {React.ReactElement} - The component showing the home page content.
 */
const Home = () => {
    return (
        <Box sx={{ maxWidth: 880, mx: 'auto', px: { xs: 2, sm: 4 }, py: { xs: 4, sm: 6 } }}>
            {/* hero */}
            <Box sx={{ textAlign: 'center' }}>
                <Box
                    sx={{
                        display: 'flex',
                        flexDirection: { xs: 'column', sm: 'row' },
                        justifyContent: 'center',
                        alignItems: 'center',
                        gap: 2,
                    }}
                >
                    <Box component='img' src='/paras.png' alt='PARAS logo' sx={LOGO_PLATE} />
                    <Box component='img' src='/parasect.png' alt='PARASECT logo' sx={LOGO_PLATE} />
                </Box>

                <Typography variant='h3' component='h1' sx={{ mt: 3, fontSize: { xs: '2rem', sm: '2.75rem' } }}>
                    PARAS &amp; PARASECT
                </Typography>
                <Typography
                    variant='h6'
                    component='p'
                    color='textSecondary'
                    sx={{ mt: 1, maxWidth: 560, mx: 'auto', fontWeight: 400 }}
                >
                    Predict the substrate specificity of adenylation domains in non-ribosomal peptide synthetases.
                </Typography>
            </Box>

            {/* the main call to action, deliberately a different shape from the rest */}
            <Card
                sx={{
                    mt: 5,
                    backgroundColor: 'primary.main',
                    color: 'primary.contrastText',
                    borderRadius: 3,
                }}
            >
                <CardActionArea
                    component={RouterLink}
                    to={primaryAction.to}
                    sx={{ p: { xs: 2.5, sm: 3 }, display: 'flex', alignItems: 'center', gap: 2 }}
                >
                    <primaryAction.Icon sx={{ fontSize: 34, flexShrink: 0 }} />
                    <Box sx={{ flex: 1, minWidth: 0 }}>
                        <Typography variant='h6' sx={{ fontWeight: 700 }}>
                            {primaryAction.label}
                        </Typography>
                        <Typography variant='body2' sx={{ opacity: 0.9 }}>
                            {primaryAction.blurb}
                        </Typography>
                    </Box>
                    <ArrowForwardIcon sx={{ flexShrink: 0, display: { xs: 'none', sm: 'block' } }} />
                </CardActionArea>
            </Card>

            {/* one grid of equal columns, so no card can size itself to its text */}
            <Box
                sx={{
                    mt: 2,
                    display: 'grid',
                    gridTemplateColumns: { xs: '1fr', sm: 'repeat(2, 1fr)' },
                    gap: 2,
                }}
            >
                {secondaryActions.map(({ to, Icon, label, blurb }) => (
                    <Card key={to} variant='outlined' sx={{ borderRadius: 3 }}>
                        {/* full-height action area: the whole card is the hit target,
                            and every card in a row ends up the same height */}
                        <CardActionArea
                            component={RouterLink}
                            to={to}
                            sx={{
                                height: '100%',
                                p: 2.5,
                                display: 'flex',
                                flexDirection: 'column',
                                alignItems: 'flex-start',
                                justifyContent: 'flex-start',
                            }}
                        >
                            <Icon color='primary' sx={{ fontSize: 26, mb: 1 }} />
                            <Typography variant='subtitle1' sx={{ fontWeight: 700 }}>
                                {label}
                            </Typography>
                            <Typography variant='body2' color='textSecondary' sx={{ mt: 0.5 }}>
                                {blurb}
                            </Typography>
                        </CardActionArea>
                    </Card>
                ))}
            </Box>

            <Divider sx={{ mt: 5, mb: 3 }} />

            <Box sx={{ display: 'flex', flexDirection: 'column', gap: 0.5, textAlign: 'center' }}>
                {links.map((link) => (
                    <Typography key={link.cta} variant='body2' color='textSecondary'>
                        {link.text}{' '}
                        <Tooltip title={link.tooltip} arrow>
                            <Link
                                href={link.href}
                                underline='hover'
                                target='_blank'
                                rel='noopener noreferrer'
                                sx={{ fontWeight: 'bold', whiteSpace: 'nowrap' }}
                            >
                                {link.cta}
                                <ExitIcon fontSize='inherit' sx={{ ml: 0.5, verticalAlign: 'text-bottom' }} />
                            </Link>
                        </Tooltip>
                    </Typography>
                ))}
            </Box>
        </Box>
    );
};

export default Home;
