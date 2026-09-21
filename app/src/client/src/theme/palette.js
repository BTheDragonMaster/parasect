/**
 * Colour tokens for the app, in one place so light and dark stay in step.
 *
 * Data colour is a separate job. Identity-by-hue has to survive red-green
 * colour blindness, which rules out picking pretty hues by hand, so CATEGORICAL
 * below is a fixed eight-slot order validated under protanopia and deuteranopia
 * simulation (Machado-Oliveira-Fernandes 2009) rather than chosen by eye:
 * worst adjacent pair dE 9.2 (target >= 8) and 19.3 for normal vision
 * (floor >= 15), in both modes. Orange leads because it is the brand hue.
 *
 * Two rules keep that guarantee intact:
 *   - slots are assigned in order and never cycled, and
 *   - a category keeps its slot for as long as it is on screen, so adding or
 *     removing one never repaints the others.
 * Past eight simultaneous categories no palette is safe. See CATEGORICAL_TIER2.
 */

/** Sampled from the logos. */
export const BRAND = {
    orange: '#F28732',   // crab shell
    red: '#B72726',      // mushroom cap
    yellow: '#F7D350',   // cap spots
    skyBlue: '#D5EDFA',  // eyes
    ink: '#1C1D1A',      // outlines
};

/** Brand orange stepped for contrast: 800 carries white text, 400 carries ink. */
const ORANGE = { 300: '#F6B67F', 400: '#F2A063', 500: BRAND.orange, 600: '#D9701F', 700: '#C2651F', 800: '#A2530F' };
const RED = { 400: '#E8615F', 500: '#C63534', 600: BRAND.red, 700: '#A02120' };

/**
 * The eight categorical slots, in fixed order, stepped per mode.
 *
 * Only the first three clear the all-pairs gate (scatter-like views, where any
 * two marks can end up side by side; the network graph is one). Beyond three
 * simultaneous categories, identity needs a second channel: the legend names
 * them, and the node labels can be turned on.
 */
export const CATEGORICAL = {
    light: ['#eb6834', '#1baf7a', '#2a78d6', '#eda100', '#e87ba4', '#008300', '#4a3aa7', '#e34948'],
    dark: ['#d95926', '#199e70', '#3987e5', '#c98500', '#d55181', '#008300', '#9085e9', '#e66767'],
};

/**
 * Slots 9-16: the same eight hues at a second lightness step.
 *
 * These are deliberately NOT validated as distinct from their tier-1 partner.
 * They differ mostly in lightness, and nothing can make sixteen hues pairwise
 * safe. They exist so a legend with more than eight entries still renders, and
 * the UI says as much when it starts using them.
 */
export const CATEGORICAL_TIER2 = {
    light: ['#9c3d18', '#0f6f4d', '#194a85', '#8f6200', '#b44a70', '#004f00', '#2d2266', '#8f2c2b'],
    dark: ['#f2a077', '#67c9a5', '#8fbcf0', '#e8bf5c', '#eda3bd', '#4db84d', '#bcb4f2', '#f0a09f'],
};

/** Nodes/rows whose category is not in the legend. Reads as "not one of the named ones". */
export const OTHER = { light: '#b8b4ac', dark: '#5c5a54' };

/** Dimmed-out marks when a highlight is active. */
export const MUTED_MARK = { light: '#dedbd4', dark: '#3a3936' };

/**
 * Signatures the user placed in the network. Every hue is already represented
 * by CATEGORICAL, and the greys by OTHER and MUTED_MARK, so a placed signature
 * is drawn in the text colour: it reads as "yours" against all of them, in
 * both modes, and never as a legend entry.
 */
export const QUERY_MARK = { light: BRAND.ink, dark: '#F2F1EE' };

const NEUTRAL_LIGHT = {
    bg: '#F7F6F3',
    paper: '#FFFFFF',
    sunken: '#EFEDE8',
    border: '#DDD9D1',
    borderStrong: '#C4BFB4',
    text: BRAND.ink,
    textSecondary: '#5A5854',
};

const NEUTRAL_DARK = {
    bg: '#191917',
    paper: '#212120',
    sunken: '#2A2A27',
    border: '#3A3935',
    borderStrong: '#4E4C47',
    text: '#F2F1EE',
    textSecondary: '#AFADA6',
};

/**
 * Per-gene banding on the results and annotation pages.
 */
export const GENE_BANDS = {
    light: [
        { surface: '#F7E6D8', rail: '#B84A15' },  // rail 4.3:1 on its band, 4.8:1 on the page
        { surface: '#E2EBF7', rail: '#1B5FA8' },  // 5.4:1 / 6.0:1
    ],
    dark: [
        { surface: '#33261E', rail: '#E8834A' },  // 5.4:1 / 6.5:1
        { surface: '#1F2A3A', rail: '#6BAEE8' },  // 6.1:1 / 7.4:1
    ],
};

export const NEUTRALS = { light: NEUTRAL_LIGHT, dark: NEUTRAL_DARK };

/** Chrome colours per mode: brand-derived, each checked for text contrast. */
export const CHROME = {
    light: {
        primary: ORANGE[800],      // white text 5.6:1
        primaryHover: '#8B460C',
        onPrimary: '#FFFFFF',
        secondary: RED[600],       // white text 6.3:1
        onSecondary: '#FFFFFF',
        accent: BRAND.orange,      // ink text 6.7:1
        onAccent: BRAND.ink,
    },
    dark: {
        primary: ORANGE[400],      // ink text 8.1:1, 8.4:1 on the dark page
        primaryHover: ORANGE[300],
        onPrimary: BRAND.ink,
        secondary: RED[400],       // 5.3:1 on the dark page
        onSecondary: BRAND.ink,
        accent: ORANGE[300],
        onAccent: BRAND.ink,
    },
};

/**
 * Signature-distance ramp: five ordinal steps for a Hamming distance.
 *
 * This is magnitude, not identity, so it is one hue stepped by lightness rather
 * than a categorical slot. And it is deliberately the brand orange rather
 * than a hue held out from CATEGORICAL, because no eight-hue set leaves one
 * spare. What keeps a distance badge from reading as "legend entry 1" is its
 * shape: it carries the number inside it, and a legend swatch never does.
 *
 * The ramp follows the number it labels: more mismatches sit further from the
 * page. Light mode runs pale -> dark, dark mode dark -> pale, so "distant" is
 * always the loudest step and "identical" the quietest. Steps are monotone in
 * OKLCH L with gaps >= 0.073, single hue (spread 1 deg), and the quiet end
 * clears 2.1:1 against its surface. Each step names the ink that clears 4.5:1
 * on it. The ramp spans the range where that flips from dark to light text.
 */
export const DISTANCE_RAMP = {
    // on paper (#FFFFFF): marks 2.07 / 2.75 / 3.63 / 5.03 / 7.52 : 1
    light: [
        { fill: '#f99f5e', text: BRAND.ink },   // text 8.18:1
        { fill: '#e9812b', text: BRAND.ink },   // 6.16:1
        { fill: '#c96f24', text: BRAND.ink },   // 4.66:1
        { fill: '#a85b1c', text: '#FFFFFF' },   // 5.03:1
        { fill: '#814512', text: '#FFFFFF' },   // 7.52:1
    ],
    // on paper (#212120): marks 2.10 / 3.07 / 4.38 / 6.19 / 9.02 : 1
    dark: [
        { fill: '#7f4412', text: '#F2F1EE' },   // text 6.79:1
        { fill: '#a3591b', text: '#F2F1EE' },   // 4.65:1
        { fill: '#c86e24', text: BRAND.ink },   // 4.61:1
        { fill: '#ef852d', text: BRAND.ink },   // 6.50:1
        { fill: '#fab281', text: BRAND.ink },   // 9.48:1
    ],
};

/** Status colours, reserved, and never reused as a categorical slot. */
export const STATUS = {
    light: { good: '#1F7A4D', warning: '#A2530F', critical: '#B3261E', info: '#1B5FA8' },
    dark: { good: '#4FBF8B', warning: '#E8A24C', critical: '#F2807A', info: '#6BAEE8' },
};

/**
 * Colour for categorical slot index, wrapping into the second tier past eight.
 *
 * @param {number} index - zero-based slot number.
 * @param {'light'|'dark'} mode - active colour mode.
 * @returns {string} - hex colour.
 */
export function categoricalColor(index, mode) {
    const tier1 = CATEGORICAL[mode];
    const tier2 = CATEGORICAL_TIER2[mode];
    if (index < tier1.length) return tier1[index];
    const rest = index - tier1.length;
    if (rest < tier2.length) return tier2[rest];
    // past sixteen the legend is far beyond anything readable; keep it stable
    // and deterministic rather than inventing more hues
    return tier2[rest % tier2.length];
}

/** How many slots carry the validated, colour-blind-safe hues. */
export const SAFE_SLOTS = CATEGORICAL.light.length;
