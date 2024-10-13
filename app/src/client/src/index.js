import React, { useState } from 'react';
import ReactDOM from 'react-dom/client';
import { BrowserRouter, Routes, Route, useNavigate } from 'react-router-dom';
import { createTheme, ThemeProvider } from '@mui/material/styles';
import { AppBar, Toolbar, IconButton, Typography, Menu, MenuItem } from '@mui/material';
import MenuIcon from '@mui/icons-material/Menu';
import HomeIcon from '@mui/icons-material/Home';
import UploadIcon from '@mui/icons-material/Upload';
import GitHubIcon from '@mui/icons-material/GitHub';

import 'react-widgets/scss/styles.scss';
import './style/main.scss';

import Toast from './components/Toast';
import Home from './pages/home';
import Submit from './pages/submit';
import Results from './pages/results';
import NotFound from './pages/not_found';

const theme = createTheme({});

const CustomToolbar = () => {
    // handle navigation
    const navigate = useNavigate();

    // version of the app
    const [version, setVersion] = useState('');

    // fetch version from server
    fetch('/api/version')
        .then((response) => response.json())
        .then((data) => setVersion(`v${data.version}`));

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
            <Toolbar>
                {/* hamburger Menu Icon */}
                <IconButton 
                    edge='start' 
                    color='inherit' 
                    aria-label='menu' 
                    onClick={handleMenuOpen}
                >
                    <MenuIcon />
                </IconButton>

                {/* menu that opens when hamburger icon is clicked */}
                <Menu
                    anchorEl={anchorEl}
                    open={open}
                    onClose={handleMenuClose}
                >
                    <MenuItem onClick={() => handleMenuItemClick('/')}>
                        <HomeIcon style={{ marginRight: '10px' }} />
                        Home
                    </MenuItem>
                    <MenuItem onClick={() => handleMenuItemClick('/submit')}>
                        <UploadIcon style={{ marginRight: '10px' }} />
                        Submit
                    </MenuItem>
                    <MenuItem onClick={() => handleExternalLinkClick('https://github.com/BTheDragonMaster/parasect/issues')}>
                        <GitHubIcon style={{ marginRight: '10px' }} />
                        Report an issue
                    </MenuItem>
                </Menu>

                {/* display name and version next to hamburger */}
                <Typography variant='h6' style={{ marginLeft: '16px' }}>
                    PARAS {version}
                </Typography>
            </Toolbar>
        </AppBar>
    );
};

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
                    path='/results/:jobId' 
                    element={<Results />}
                />
                <Route 
                    path='*' 
                    element={<NotFound />}
                />
            </Routes>
        </div>
    );
};

function App () {
    return (
        <ThemeProvider theme={theme}>
            <BrowserRouter>
                <CustomToolbar />
                <AppRoutes />
                <Toast />
            </BrowserRouter>
        </ThemeProvider>
    );
};

const root = ReactDOM.createRoot(document.getElementById('root'));

root.render(<App />);