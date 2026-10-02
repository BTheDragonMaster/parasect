/**
 * Helpers for the signature comparison page: colouring residues by their
 * similarity to the reference row, reordering rows, and the URL other pages use
 * to open the comparison prefilled.
 */

/** Which signature a row is compared on. */
export const SIGNATURE_MODES = {
    extended: { label: 'Extended signature (34)', field: 'extended_signature' },
    signature: { label: 'Signature (10)', field: 'signature' },
};

/**
 * Label columns that can be shown next to the sequences, in menu order.
 * `get` reads the value from a row; custom (pasted) rows have none of the
 * database fields and show an empty cell.
 */
export const LABEL_COLUMNS = {
    name: { label: 'Domain', get: (r) => r.name },
    protein: { label: 'Protein', get: (r) => r.protein },
    substrates: { label: 'Substrate', get: (r) => (r.substrates || []).join(', ') },
    species: { label: 'Species', get: (r) => r.taxonomy?.species },
    strain: { label: 'Strain', get: (r) => r.taxonomy?.strain },
    genus: { label: 'Genus', get: (r) => r.taxonomy?.genus },
    family: { label: 'Family', get: (r) => r.taxonomy?.family },
    order: { label: 'Order', get: (r) => r.taxonomy?.order },
    cls: { label: 'Class', get: (r) => r.taxonomy?.cls },
    phylum: { label: 'Phylum', get: (r) => r.taxonomy?.phylum },
    kingdom: { label: 'Kingdom', get: (r) => r.taxonomy?.kingdom },
    domain: { label: 'Domain of life', get: (r) => r.taxonomy?.domain },
};

export const DEFAULT_LABEL_COLUMNS = ['name', 'substrates', 'genus'];

const hexToRgb = (hex) => {
    const n = parseInt(hex.replace('#', ''), 16);
    return [(n >> 16) & 255, (n >> 8) & 255, n & 255];
};

const rgbToHex = (rgb) => `#${rgb.map((c) => Math.round(c).toString(16).padStart(2, '0')).join('')}`;

/** WCAG relative luminance of an sRGB colour. */
const luminance = (rgb) => {
    const [r, g, b] = rgb.map((c) => {
        const v = c / 255;
        return v <= 0.03928 ? v / 12.92 : ((v + 0.055) / 1.055) ** 2.4;
    });
    return 0.2126 * r + 0.7152 * g + 0.0722 * b;
};

const contrast = (a, b) => {
    const [hi, lo] = [luminance(a), luminance(b)].sort((x, y) => y - x);
    return (hi + 0.05) / (lo + 0.05);
};

/**
 * Fill and text colour for a residue with the given similarity to the reference.
 *
 * @param {number} similarity - in [0, 1]; 1 is identical.
 * @param {{low: string, high: string, inks: string[]}} ramp - SIMILARITY_RAMP for the mode.
 * @returns {{fill: string, text: string}} - cell colours.
 */
export function similarityColors(similarity, ramp) {
    const t = Math.min(Math.max(similarity, 0), 1);
    const low = hexToRgb(ramp.low);
    const high = hexToRgb(ramp.high);
    const fill = low.map((c, i) => c + (high[i] - c) * t);
    const inks = ramp.inks.map(hexToRgb);
    const text = contrast(fill, inks[0]) >= contrast(fill, inks[1]) ? ramp.inks[0] : ramp.inks[1];
    return { fill: rgbToHex(fill), text };
}

/**
 * Look up how alike two residues are.
 *
 * @param {{alphabet: string, matrix: number[][]}|null} table - from /api/compare/aa_similarity.
 * @param {string} a - residue.
 * @param {string} b - residue.
 * @returns {number|null} - similarity in [0, 1], or null for a residue the table doesn't know.
 */
export function residueSimilarity(table, a, b) {
    if (!table || !a || !b) return null;
    const i = table.alphabet.indexOf(a.toUpperCase());
    const j = table.alphabet.indexOf(b.toUpperCase());
    if (i < 0 || j < 0) return null;
    return table.matrix[i][j];
}

/**
 * Move one item of a list to another index.
 *
 * @param {Array} list - the list.
 * @param {number} from - index of the item to move.
 * @param {number} to - index it should end up at.
 * @returns {Array} - a new list.
 */
export function moveItem(list, from, to) {
    if (from === to || from < 0 || from >= list.length) return list;
    const next = [...list];
    const [item] = next.splice(from, 1);
    next.splice(Math.max(0, Math.min(to, next.length)), 0, item);
    return next;
}

/**
 * Link to the comparison page, prefilled.
 *
 * Custom (non-database) signatures come first, so a placed query is the
 * reference its neighbours are compared against; database domains follow in
 * the order given.
 *
 * @param {{ids?: number[], custom?: Array<{name: string, signature: string}>}} entries - what to show.
 * @returns {string} - a path with query string.
 */
export function compareUrl({ ids = [], custom = [] }) {
    const params = new URLSearchParams();
    custom.forEach(({ name, signature }) => params.append('custom', `${name}:${signature}`));
    if (ids.length) params.set('ids', ids.join(','));
    return `/compare?${params}`;
}

/**
 * Read the prefill back out of the page's query string.
 *
 * @param {URLSearchParams} params - the page's search params.
 * @returns {{ids: number[], custom: Array<{name: string, signature: string}>}} - what to load.
 */
export function parseCompareParams(params) {
    const ids = (params.get('ids') || '')
        .split(',')
        .map((v) => parseInt(v, 10))
        .filter((v) => Number.isInteger(v));
    const custom = params.getAll('custom').map((value) => {
        // the name may itself contain a colon; the signature never does
        const at = value.lastIndexOf(':');
        return at < 0
            ? { name: 'Custom signature', signature: value }
            : { name: value.slice(0, at) || 'Custom signature', signature: value.slice(at + 1) };
    }).filter((c) => c.signature);
    return { ids, custom };
}

/**
 * Hamming distance between two position-matched signatures.
 *
 * Signatures of one kind all have the same length; should two differ anyway,
 * each missing position counts as a mismatch.
 *
 * @param {string} a - signature.
 * @param {string} b - signature.
 * @returns {number|null} - mismatching positions, or null if either is missing.
 */
export function hammingDistance(a, b) {
    if (!a || !b) return null;
    let distance = Math.abs(a.length - b.length);
    for (let i = 0; i < Math.min(a.length, b.length); i += 1) {
        if (a[i] !== b[i]) distance += 1;
    }
    return distance;
}

/**
 * Order rows by Hamming distance to the first (reference) row, which stays on top.
 *
 * The sort is stable, so rows at the same distance keep their current order,
 * and rows without this kind of signature go last either way.
 *
 * @param {Array<object>} rows - rows, reference first.
 * @param {string} field - 'signature' or 'extended_signature'.
 * @param {'asc'|'desc'} direction - closest first, or furthest first.
 * @returns {Array<object>} - a new list.
 */
export function sortByDistance(rows, field, direction) {
    if (rows.length < 3) return rows;
    const [reference, ...rest] = rows;
    const sign = direction === 'desc' ? -1 : 1;
    const keyed = rest.map((row) => ({ row, d: hammingDistance(reference[field], row[field]) }));
    keyed.sort((x, y) => {
        if (x.d === null || y.d === null) return (x.d === null) - (y.d === null);
        return sign * (x.d - y.d);
    });
    return [reference, ...keyed.map((k) => k.row)];
}
