import React, { createContext, useCallback, useContext, useEffect, useMemo, useState } from 'react';
import { ThemeProvider } from '@mui/material/styles';
import CssBaseline from '@mui/material/CssBaseline';

import { createAppTheme } from './index';

const STORAGE_KEY = 'paras-color-mode';
const MODES = ['light', 'dark', 'system'];

const ColorModeContext = createContext({
    /** what the user picked: 'light' | 'dark' | 'system' */
    preference: 'system',
    /** what that resolves to right now: 'light' | 'dark' */
    mode: 'light',
    setPreference: () => {},
});

/** @returns {{preference: string, mode: string, setPreference: Function}} - the colour-mode controls. */
export const useColorMode = () => useContext(ColorModeContext);

const prefersDark = () => typeof window !== 'undefined'
    && window.matchMedia?.('(prefers-color-scheme: dark)').matches;

function readStoredPreference() {
    try {
        const stored = window.localStorage.getItem(STORAGE_KEY);
        return MODES.includes(stored) ? stored : 'system';
    } catch (err) {
        // private browsing or blocked storage: fall back to following the OS
        return 'system';
    }
}

/**
 * Supplies the theme and the light/dark/system control to the whole app.
 *
 * 'system' is the default and stays live: the OS listener below re-resolves the
 * mode when the user flips their machine to dark at dusk, without a reload.
 *
 * @param {object} props - component props.
 * @param {React.ReactNode} props.children - the app tree.
 * @returns {React.ReactElement} - the themed tree.
 */
export function ColorModeProvider({ children }) {
    const [preference, setPreferenceState] = useState(readStoredPreference);
    const [systemIsDark, setSystemIsDark] = useState(prefersDark);

    useEffect(() => {
        const query = window.matchMedia?.('(prefers-color-scheme: dark)');
        if (!query) return undefined;
        const onChange = (event) => setSystemIsDark(event.matches);
        query.addEventListener('change', onChange);
        return () => query.removeEventListener('change', onChange);
    }, []);

    const setPreference = useCallback((next) => {
        if (!MODES.includes(next)) return;
        setPreferenceState(next);
        try {
            window.localStorage.setItem(STORAGE_KEY, next);
        } catch (err) {
            // the preference still applies for this session; nothing to recover
        }
    }, []);

    const mode = preference === 'system' ? (systemIsDark ? 'dark' : 'light') : preference;
    const theme = useMemo(() => createAppTheme(mode), [mode]);

    // let plain CSS (the toast styles, scrollbars) react to the mode too
    useEffect(() => {
        document.documentElement.setAttribute('data-theme', mode);
        document.documentElement.style.colorScheme = mode;
    }, [mode]);

    const value = useMemo(() => ({ preference, mode, setPreference }), [preference, mode, setPreference]);

    return (
        <ColorModeContext.Provider value={value}>
            <ThemeProvider theme={theme}>
                <CssBaseline />
                {children}
            </ThemeProvider>
        </ColorModeContext.Provider>
    );
}

export default ColorModeContext;
