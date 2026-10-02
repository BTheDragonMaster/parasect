/**
 * Text helpers shared by the SVG figure exports.
 */

export const FONT = "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif";
export const MONO_FONT = "ui-monospace, SFMono-Regular, Menlo, Consolas, monospace";

export const escapeXml = (value) => String(value)
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;');

/** Round to 2dp; figure coordinates need no more. */
export const r2 = (n) => Math.round(n * 100) / 100;

let measureContext = null;

/**
 * Width of a string in the figure font, measured on a canvas; falls back to an
 * estimate where there is no canvas (tests).
 */
export const textWidth = (text, size, weight = 400, font = FONT) => {
    if (measureContext === null && typeof document !== 'undefined') {
        measureContext = document.createElement('canvas').getContext('2d') || undefined;
    }
    if (!measureContext) return String(text).length * size * 0.56;
    measureContext.font = `${weight} ${size}px ${font}`;
    return measureContext.measureText(String(text)).width;
};

/** Cut `text` down, with an ellipsis, until it fits in `maxWidth`. */
export const fitText = (text, maxWidth, size, weight = 400) => {
    const full = String(text);
    if (textWidth(full, size, weight) <= maxWidth) return full;
    let lo = 0;
    let hi = full.length;
    while (lo < hi) {
        const mid = Math.ceil((lo + hi) / 2);
        if (textWidth(`${full.slice(0, mid)}…`, size, weight) <= maxWidth) lo = mid;
        else hi = mid - 1;
    }
    return lo > 0 ? `${full.slice(0, lo)}…` : '…';
};
