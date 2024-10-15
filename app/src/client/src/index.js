import React, { useState } from 'react';
import ReactDOM from 'react-dom/client';
import { BrowserRouter, Routes, Route, useNavigate } from 'react-router-dom';
import { createTheme, ThemeProvider } from '@mui/material/styles';
import { AppBar, Toolbar, IconButton, Typography, Menu, MenuItem } from '@mui/material';
import { MdMenu } from 'react-icons/md';
import HomeIcon from '@mui/icons-material/Home';
import UploadIcon from '@mui/icons-material/Upload';
import GitHubIcon from '@mui/icons-material/GitHub';

import './style/main.css';

import Toast from './components/Toast';
import Home from './pages/home';
import Submit from './pages/submit';
import Results from './pages/results';
import NotFound from './pages/not_found';

/**
 * Custom theme for the app.
 * 
 * @returns {Theme} - The custom theme for the app.
 */
const theme = createTheme({
    palette: {
        primary: {
            main: '#3d7dca',
        },
        secondary: {
            main: '#ffcb05',
        },
        white: {
            main: '#ffffff',
        },
        black: {
            main: '#000000',
        },
        gray: {
            main: '#f5f5f5',
        },
    },
    typography: {
        fontFamily: [
            'Arial',
            'Roboto',
            'sans-serif',
        ].join(','),
    },
});

/**
 * Custom toolbar for the app.
 * 
 * @returns {React.ReactElement} - The custom toolbar for the app.
 */
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
        <AppBar position='static' sx={{ backgroundColor: 'primary.main' }}>
            <Toolbar>
                {/* hamburger Menu Icon */}
                <IconButton onClick={handleMenuOpen} sx={{ mr: 2}}>
                    <MdMenu fill='white' />
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
                    <MenuItem onClick={() => handleExternalLinkClick('https://github.com/BTheDragonMaster/parasect/issues')}>
                        <GitHubIcon sx={{ marginRight: '10px' }} />
                        Report an issue
                    </MenuItem>
                </Menu>

                {/* display name and version next to hamburger */}
                <Typography 
                    variant='h6' 
                    sx={{ marginLeft: '16px' }}
                >
                    <Typography sx={{ color: 'white.main' }}>
                        PARAS {version}
                    </Typography>
                </Typography>
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

/**
 * Main app component.
 * 
 * @returns {React.ReactElement} - The main app component.
 */
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
