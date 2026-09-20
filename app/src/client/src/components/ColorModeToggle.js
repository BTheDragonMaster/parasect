import React, { useState } from 'react';
import { IconButton, Menu, MenuItem, ListItemIcon, ListItemText, Tooltip } from '@mui/material';
import LightModeIcon from '@mui/icons-material/LightMode';
import DarkModeIcon from '@mui/icons-material/DarkMode';
import SettingsBrightnessIcon from '@mui/icons-material/SettingsBrightness';
import CheckIcon from '@mui/icons-material/Check';

import { useColorMode } from '../theme/ColorModeContext';

const OPTIONS = [
    { value: 'light', label: 'Light', Icon: LightModeIcon },
    { value: 'dark', label: 'Dark', Icon: DarkModeIcon },
    { value: 'system', label: 'System', Icon: SettingsBrightnessIcon },
];

/**
 * App-bar control for the light / dark / system colour mode.
 *
 * @returns {React.ReactElement} - the toggle.
 */
const ColorModeToggle = () => {
    const { preference, mode, setPreference } = useColorMode();
    const [anchorEl, setAnchorEl] = useState(null);

    const active = OPTIONS.find((o) => o.value === preference) || OPTIONS[2];
    // the button shows what you are looking at, not what you picked, so
    // 'system' still reads as a sun or a moon
    const CurrentIcon = mode === 'dark' ? DarkModeIcon : LightModeIcon;

    return (
        <>
            <Tooltip title={`Appearance: ${active.label.toLowerCase()}`}>
                <IconButton
                    onClick={(e) => setAnchorEl(e.currentTarget)}
                    size='small'
                    aria-label='Change appearance'
                    sx={{ color: 'inherit' }}
                >
                    <CurrentIcon fontSize='small' />
                </IconButton>
            </Tooltip>
            <Menu anchorEl={anchorEl} open={Boolean(anchorEl)} onClose={() => setAnchorEl(null)}>
                {OPTIONS.map(({ value, label, Icon }) => (
                    <MenuItem
                        key={value}
                        selected={preference === value}
                        onClick={() => { setPreference(value); setAnchorEl(null); }}
                    >
                        <ListItemIcon><Icon fontSize='small' /></ListItemIcon>
                        <ListItemText>{label}</ListItemText>
                        {preference === value && <CheckIcon fontSize='small' sx={{ ml: 2, opacity: 0.7 }} />}
                    </MenuItem>
                ))}
            </Menu>
        </>
    );
};

export default ColorModeToggle;
