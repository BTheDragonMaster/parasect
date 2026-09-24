import React, { useEffect, useState } from 'react';
import ReactDOM from 'react-dom/client';
import { BrowserRouter, Routes, Route, useNavigate } from 'react-router-dom';
import { AppBar, Toolbar, Box, IconButton, Typography, Menu, MenuItem } from '@mui/material';
import { MdMenu } from 'react-icons/md';
import HomeIcon from '@mui/icons-material/Home';
import UploadIcon from '@mui/icons-material/Upload';
import GitHubIcon from '@mui/icons-material/GitHub';
import RetrieveIcon from '@mui/icons-material/GetApp';
import DatasetIcon from '@mui/icons-material/Dataset'
import QueryStatsIcon from '@mui/icons-material/QueryStats';
import HubIcon from '@mui/icons-material/Hub';
import StorageIcon from '@mui/icons-material/Storage';

import './style/main.css';

import { GITHUB_ISSUES_URL } from './links';
import { ColorModeProvider } from './theme/ColorModeContext';
import ColorModeToggle from './components/ColorModeToggle';
import Toast from './components/Toast';
import Home from './pages/home';
import Retrieve from './pages/retrieve';
import Submit from './pages/submit';
import Results from './pages/results';
import NotFound from './pages/not_found';
import DataAnnotation from './pages/data_annotation'
import AnnotationEditor from './pages/annotation_editor'
import QueryDatabase from './pages/query_database'
import NetworkGraph from './pages/network'
import Dataset from './pages/dataset'

/**
 * Custom toolbar for the app.
 * 
 * @returns {React.ReactElement} - The custom toolbar for the app.
 */
const CustomToolbar = () => {
    // handle navigation
    const navigate = useNavigate();

    // version of the app
    const [version, setVersion] = useState('UNKNOWN');

    // fetch version from server
    useEffect(() => {
    fetch('/api/version')
        .then((response) => {
            if (!response.ok) {
                throw new Error(`HTTP error! status: ${response.status}`);
            }
            return response.json();
        })
        .then((data) => {
            setVersion(`v${data.version}`);
        })
        .catch((error) => {
            console.error('Failed to fetch version:', error);
            setVersion('v?');  // Fallback version if fetch fails
        });
}, []);

    // state to handle menu
    const [anchorEl, setAnchorEl] = useState(null);
    const open = Boolean(anchorEl);

    // function to handle opening menu
    const handleMenuOpen = (event) => {
        setAnchorEl(event.currentTarget);
    };

    // function to handle closing menu
    const handleMenuClose = () => {
        setAnchorEl(null);
    };

    // handle meny item click
    const handleMenuItemClick = (path) => {
        navigate(path);
        handleMenuClose();
    };

    // function to open external link in a new tab
    const handleExternalLinkClick = (url) => {
        window.open(url, '_blank');  // open the link in a new tab
        handleMenuClose();  // close the menu after opening
    };

    return (
        <AppBar position='static'>
            <Toolbar sx={{ gap: 1 }}>
                {/* hamburger Menu Icon */}
                <IconButton onClick={handleMenuOpen} edge='start' sx={{ mr: 1, color: 'inherit' }}>
                    <MdMenu fill='currentColor' size={22} />
                </IconButton>

                {/* menu that opens when hamburger icon is clicked */}
                <Menu
                    anchorEl={anchorEl}
                    open={open}
                    onClose={handleMenuClose}
                >
                    <MenuItem onClick={() => handleMenuItemClick('/')}>
                        <HomeIcon sx={{ marginRight: '10px' }} />
                        Home
                    </MenuItem>
                    <MenuItem onClick={() => handleMenuItemClick('/submit')}>
                        <UploadIcon sx={{ marginRight: '10px' }} />
                        Submit
                    </MenuItem>
                    <MenuItem onClick={() => handleMenuItemClick('/retrieve')}>
                        <RetrieveIcon sx={{ marginRight: '10px' }} />
                        Retrieve
                    </MenuItem>
                    <MenuItem onClick={() => handleMenuItemClick('/data_annotation')}>
                        <DatasetIcon sx={{ marginRight: '10px' }} />
                        Data annotation
                    </MenuItem>
                    <MenuItem onClick={() => handleMenuItemClick('/query_database')}>
                        <QueryStatsIcon sx={{ marginRight: '10px' }} />
                        Query database
                    </MenuItem>
                    <MenuItem onClick={() => handleMenuItemClick('/network')}>
                        <HubIcon sx={{ marginRight: '10px' }} />
                        Network
                    </MenuItem>
                    <MenuItem onClick={() => handleMenuItemClick('/dataset')}>
                        <StorageIcon sx={{ marginRight: '10px' }} />
                        Download dataset
                    </MenuItem>
                    <MenuItem onClick={() => handleExternalLinkClick(GITHUB_ISSUES_URL)}>
                        <GitHubIcon sx={{ marginRight: '10px' }} />
                        Report an issue
                    </MenuItem>
                </Menu>

                {/* display name and version next to hamburger */}
                <Box sx={{ minWidth: 0, overflow: 'hidden' }}>
                    <Typography
                        noWrap
                        sx={{
                            color: 'inherit',
                            fontWeight: 600,
                            fontSize: { xs: '0.95rem', sm: '1.15rem' },
                        }}
                    >
                        PARAS {version}
                    </Typography>
                    <Typography
                        noWrap
                        sx={{
                            color: 'inherit',
                            opacity: 0.8,
                            fontSize: '0.7rem',
                            display: { xs: 'none', sm: 'block' },
                        }}
                    >
                        web app v{process.env.REACT_APP_VERSION ? process.env.REACT_APP_VERSION : 'UNKNOWN'}
                    </Typography>
                </Box>

                <Box sx={{ flexGrow: 1 }} />
                <ColorModeToggle />
            </Toolbar>
        </AppBar>
    );
};

/**
 * App routes for the app.
 * 
 * @returns {React.ReactElement} - The app routes for the app.
 */
function AppRoutes () {
    return (
        <div>
            <Routes>
                <Route 
                    path='/' 
                    element={<Home />}
                />
                <Route 
                    path='/submit' 
                    element={<Submit />}
                />
                <Route 
                    path='/retrieve' 
                    element={<Retrieve />}
                />
                <Route 
                    path='/results/:jobId' 
                    element={<Results />}
                />
                <Route
                    path='/annotation_editor/:jobId'
                    element={<AnnotationEditor />}
                />
                <Route
                    path='/data_annotation'
                    element={<DataAnnotation />}
                />
                <Route
                    path='/query_database'
                    element={<QueryDatabase />}
                />
                <Route
                    path='/network'
                    element={<NetworkGraph />}
                />
                <Route
                    path='/dataset'
                    element={<Dataset />}
                />
                <Route
                    path='*'
                    element={<NotFound />}
                />
            </Routes>
        </div>
    );
};

/**
 * Main app component.
 * 
 * @returns {React.ReactElement} - The main app component.
 */
function App () {
    return (
        <ColorModeProvider>
            <BrowserRouter>
                <CustomToolbar />
                <AppRoutes />
                <Toast />
            </BrowserRouter>
        </ColorModeProvider>
    );
};

const root = ReactDOM.createRoot(document.getElementById('root'));

root.render(<App />);
