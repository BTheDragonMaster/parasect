import { createTheme } from '@mui/material/styles';

import { BRAND, CHROME, GENE_BANDS, NEUTRALS, STATUS, categoricalColor, SAFE_SLOTS } from './palette';

export { BRAND, CHROME, GENE_BANDS, NEUTRALS, STATUS, categoricalColor, SAFE_SLOTS };
export { CATEGORICAL, CATEGORICAL_TIER2, DISTANCE_RAMP, MUTED_MARK, OTHER, QUERY_MARK, SIMILARITY_RAMP } from './palette';

const FONT_STACK = [
    '-apple-system',
    'BlinkMacSystemFont',
    '"Segoe UI"',
    'Roboto',
    'Helvetica',
    'Arial',
    'sans-serif',
].join(',');

/**
 * Build the MUI theme for a colour mode.
 *
 * Dark is a selected set of steps rather than an inversion of light: the brand
 * orange that carries white text at 5.6:1 on paper is far too dark to sit on a
 * dark page, so each mode names its own step (see CHROME in ./palette).
 *
 * @param {'light'|'dark'} mode - resolved colour mode, never 'system'.
 * @returns {import('@mui/material/styles').Theme} - the theme for that mode.
 */
export function createAppTheme(mode) {
    const n = NEUTRALS[mode];
    const c = CHROME[mode];
    const s = STATUS[mode];
    const isDark = mode === 'dark';

    return createTheme({
        palette: {
            mode,
            primary: { main: c.primary, dark: c.primaryHover, contrastText: c.onPrimary },
            secondary: { main: c.secondary, contrastText: c.onSecondary },
            // the logo orange, for identity accents rather than actions
            accent: { main: c.accent, contrastText: c.onAccent },
            success: { main: s.good },
            warning: { main: s.warning },
            error: { main: s.critical },
            info: { main: s.info },
            background: { default: n.bg, paper: n.paper },
            text: { primary: n.text, secondary: n.textSecondary },
            divider: n.border,
            // kept so the older pages that say `white.main` / `gray.main` keep
            // working; they mean "a surface", which is mode-dependent now
            white: { main: n.paper },
            black: { main: n.text },
            gray: { main: n.sunken },
            // app-specific tokens the components read directly
            surface: { sunken: n.sunken, border: n.border, borderStrong: n.borderStrong },
            geneBands: GENE_BANDS[mode],
        },
        shape: { borderRadius: 10 },
        typography: {
            fontFamily: FONT_STACK,
            h3: { fontWeight: 700 },
            h4: { fontWeight: 700 },
            h5: { fontWeight: 600 },
            h6: { fontWeight: 600 },
            button: { fontWeight: 600, textTransform: 'none' },
        },
        components: {
            MuiCssBaseline: {
                styleOverrides: {
                    // tells the browser to render form controls and scrollbars
                    // for this mode, which CssBaseline alone does not do
                    ':root': { colorScheme: mode },
                    body: { backgroundColor: n.bg, transition: 'background-color 160ms ease' },
                },
            },
            MuiButton: {
                defaultProps: { disableElevation: true },
                styleOverrides: { root: { borderRadius: 8 } },
            },
            MuiAppBar: {
                styleOverrides: {
                    root: {
                        backgroundColor: isDark ? n.paper : c.primary,
                        color: isDark ? n.text : c.onPrimary,
                        boxShadow: isDark
                            ? `0 1px 0 ${n.border}`
                            : '0 1px 3px rgba(28, 29, 26, 0.2)',
                    },
                },
            },
            MuiPaper: {
                styleOverrides: {
                    root: { backgroundImage: 'none' },
                    outlined: { borderColor: n.border },
                },
            },
            MuiDivider: { styleOverrides: { root: { borderColor: n.border } } },
            MuiTooltip: {
                styleOverrides: {
                    tooltip: {
                        backgroundColor: isDark ? n.sunken : BRAND.ink,
                        color: isDark ? n.text : '#FFFFFF',
                        fontSize: '0.75rem',
                    },
                },
            },
        },
    });
}
