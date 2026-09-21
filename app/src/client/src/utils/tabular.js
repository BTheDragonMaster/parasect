/**
 * Delimited-text export shared by the pages that offer a table download.
 */

import { downloadBlob } from './zip';

/**
 * One cell as delimited text.
 *
 * Lists are joined with ';', never ',': substrate names such as
 * "2,3-diaminopropionic acid" carry commas of their own. In TSV a stray tab or
 * line break would shift every column after it, so those become spaces; CSV
 * quotes the cell instead.
 *
 * @param {*} value - raw cell value.
 * @param {string} delim - column delimiter.
 * @returns {string} - the cell, safe to join with `delim`.
 */
function formatCell(value, delim) {
    const text = Array.isArray(value) ? value.join(';') : String(value ?? '');
    if (delim === '\t') return text.replace(/[\t\r\n]+/g, ' ');
    const flat = text.replace(/[\r\n]+/g, ' ');
    return flat.includes(delim) || flat.includes('"') ? `"${flat.replace(/"/g, '""')}"` : flat;
}

/**
 * Rows of objects as delimited text with a header line.
 *
 * @param {Array<{field: string}>} cols - columns, in order; `field` is both the key and the header.
 * @param {Array<object>} data - one object per row.
 * @param {string} delim - column delimiter, '\t' or ','.
 * @returns {string} - the table.
 */
export function makeDelimited(cols, data, delim) {
    const header = cols.map((c) => formatCell(c.field, delim)).join(delim);
    const lines = data.map((row) => cols.map((c) => formatCell(row[c.field], delim)).join(delim));
    return [header, ...lines].join('\n');
}

/**
 * Hand a string to the browser as a file download.
 *
 * @param {string} content - file contents.
 * @param {string} filename - suggested name.
 * @param {string} mime - media type.
 * @returns {void}
 */
export function downloadFile(content, filename, mime) {
    downloadBlob(new Blob([content], { type: mime }), filename);
}

export const TSV_MIME = 'text/tab-separated-values';
