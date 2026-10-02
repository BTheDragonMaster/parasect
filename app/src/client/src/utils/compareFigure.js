/**
 * Export the signature comparison as a figure: the rows in their on-screen
 * order, with the chosen label columns, each residue shaded as on the page,
 * and the distance/identity columns. Plain SVG, so it stays editable; the PNG
 * is rasterised from it (svgToPng in ./networkExport).
 */

import { FONT, MONO_FONT, escapeXml, fitText, r2, textWidth } from './svgText';

const PADDING = 24;
const HEADER_HEIGHT = 62;
const COLUMN_HEAD_HEIGHT = 22;
const ROW_HEIGHT = 22;
const CELL_WIDTH = 18;
const LABEL_SIZE = 11;
const MAX_LABEL_WIDTH = 240;
const COLUMN_GAP = 14;
const LEGEND_HEIGHT = 40;

/**
 * Render the comparison to a standalone SVG document.
 *
 * @param {object} options - export options.
 * @param {Array<{label: string, values: string[]}>} options.columns - label columns, values per row.
 * @param {Array<{cells: Array<{residue: string, fill?: string, text?: string, reference?: boolean}>,
 *     distance: (number|null), identical: number, compared: number}>} options.rows - rows in order,
 *     the reference first.
 * @param {number} options.length - the longest signature, for the position ruler.
 * @param {{title: string, subtitle: string}} options.caption - figure heading.
 * @param {{background: string, text: string, textSecondary: string, border: string, strongBorder: string,
 *     referenceFill: string, rampLow: string, rampHigh: string}} options.colors - colours.
 * @returns {string} - the SVG source.
 */
export function comparisonToSvg({ columns, rows, length, caption, colors }) {
    // each label column is as wide as its widest value, within a cap
    const labelColumns = columns.map((column) => {
        const values = column.values.map((v) => fitText(v || '', MAX_LABEL_WIDTH, LABEL_SIZE));
        const width = Math.ceil(Math.max(
            textWidth(column.label, LABEL_SIZE, 700),
            ...values.map((v, i) => textWidth(v, LABEL_SIZE, i === 0 ? 700 : 400)),
        ));
        return { label: column.label, values, width };
    });

    let x = PADDING;
    labelColumns.forEach((column) => {
        column.x = x;
        x += column.width + COLUMN_GAP;
    });
    const sequenceX = x;
    const sequenceWidth = length * CELL_WIDTH;
    const distanceX = sequenceX + sequenceWidth + COLUMN_GAP;
    const distanceWidth = Math.ceil(textWidth('Distance', LABEL_SIZE, 700));
    const identicalX = distanceX + distanceWidth + COLUMN_GAP;
    const identicalWidth = Math.ceil(Math.max(textWidth('Identical', LABEL_SIZE, 700), textWidth('00/00', LABEL_SIZE)));
    const tableRight = identicalX + identicalWidth;

    const titleWidth = Math.max(textWidth(caption.title, 19, 700), textWidth(caption.subtitle, 12));
    const totalWidth = Math.ceil(Math.max(tableRight + PADDING, PADDING * 2 + titleWidth, PADDING * 2 + 320));
    const tableTop = HEADER_HEIGHT + COLUMN_HEAD_HEIGHT;
    const tableBottom = tableTop + rows.length * ROW_HEIGHT;
    const totalHeight = tableBottom + LEGEND_HEIGHT;

    const parts = [];
    parts.push('<?xml version="1.0" encoding="UTF-8"?>');
    parts.push(`<svg xmlns="http://www.w3.org/2000/svg" width="${totalWidth}" height="${totalHeight}" `
        + `viewBox="0 0 ${totalWidth} ${totalHeight}" font-family="${FONT}">`);
    parts.push(`<defs><linearGradient id="similarity-ramp"><stop offset="0" stop-color="${colors.rampLow}"/>`
        + `<stop offset="1" stop-color="${colors.rampHigh}"/></linearGradient></defs>`);
    parts.push(`<rect width="${totalWidth}" height="${totalHeight}" fill="${colors.background}"/>`);

    parts.push(`<text x="${PADDING}" y="30" font-size="19" font-weight="700" fill="${colors.text}">`
        + `${escapeXml(caption.title)}</text>`);
    parts.push(`<text x="${PADDING}" y="50" font-size="12" fill="${colors.textSecondary}">`
        + `${escapeXml(caption.subtitle)}</text>`);

    // column heads, and a position ruler over the residues
    const headY = tableTop - 7;
    parts.push(`<g font-size="${LABEL_SIZE}" font-weight="700" fill="${colors.text}">`);
    labelColumns.forEach((column) => {
        parts.push(`<text x="${column.x}" y="${headY}">${escapeXml(column.label)}</text>`);
    });
    parts.push(`<text x="${distanceX + distanceWidth}" y="${headY}" text-anchor="end">Distance</text>`);
    parts.push(`<text x="${identicalX}" y="${headY}">Identical</text>`);
    parts.push('</g>');
    parts.push(`<g font-size="9" fill="${colors.textSecondary}" text-anchor="middle">`);
    for (let i = 0; i < length; i += 1) {
        if (i === 0 || (i + 1) % 5 === 0) {
            parts.push(`<text x="${r2(sequenceX + i * CELL_WIDTH + CELL_WIDTH / 2)}" y="${headY}">${i + 1}</text>`);
        }
    }
    parts.push('</g>');

    rows.forEach((row, index) => {
        const y = tableTop + index * ROW_HEIGHT;
        const mid = r2(y + ROW_HEIGHT / 2);
        const isReference = index === 0;
        parts.push('<g>');

        labelColumns.forEach((column) => {
            const value = column.values[index];
            parts.push(`<text x="${column.x}" y="${mid}" dominant-baseline="central" font-size="${LABEL_SIZE}" `
                + `font-weight="${isReference ? 700 : 400}" fill="${value ? colors.text : colors.textSecondary}">`
                + `${escapeXml(value || '–')}</text>`);
        });

        if (row.cells.length) {
            row.cells.forEach((cell, position) => {
                const cx = sequenceX + position * CELL_WIDTH;
                const fill = cell.reference ? colors.referenceFill : (cell.fill || colors.background);
                const ink = cell.reference ? colors.text : (cell.text || colors.text);
                parts.push(`<rect x="${cx}" y="${y + 1}" width="${CELL_WIDTH}" height="${ROW_HEIGHT - 2}" fill="${fill}"/>`);
                parts.push(`<text x="${r2(cx + CELL_WIDTH / 2)}" y="${mid}" text-anchor="middle" dominant-baseline="central" `
                    + `font-family="${MONO_FONT}" font-size="12" font-weight="${cell.reference ? 700 : 500}" fill="${ink}">`
                    + `${escapeXml(cell.residue)}</text>`);
            });
        } else {
            parts.push(`<text x="${sequenceX + 4}" y="${mid}" dominant-baseline="central" font-size="10" `
                + `fill="${colors.textSecondary}">not available</text>`);
        }

        if (!isReference && row.distance !== null) {
            parts.push(`<text x="${distanceX + distanceWidth}" y="${mid}" dominant-baseline="central" text-anchor="end" `
                + `font-size="${LABEL_SIZE}" fill="${colors.text}">${row.distance}</text>`);
        }
        if (!isReference && row.compared) {
            parts.push(`<text x="${identicalX}" y="${mid}" dominant-baseline="central" font-size="${LABEL_SIZE}" `
                + `fill="${colors.textSecondary}">${row.identical}/${row.compared}</text>`);
        }

        // the reference is set off by a heavier rule, like on screen
        parts.push(`<line x1="${PADDING}" y1="${y + ROW_HEIGHT}" x2="${tableRight}" y2="${y + ROW_HEIGHT}" `
            + `stroke="${isReference ? colors.strongBorder : colors.border}" stroke-width="${isReference ? 1.5 : 0.75}"/>`);
        parts.push('</g>');
    });

    // legend: the shading ramp
    const legendY = tableBottom + 22;
    const dissimilarWidth = textWidth('Dissimilar', 10);
    parts.push(`<g font-size="10" fill="${colors.textSecondary}">`);
    parts.push(`<text x="${PADDING}" y="${legendY}" dominant-baseline="central">Dissimilar</text>`);
    parts.push(`<rect x="${r2(PADDING + dissimilarWidth + 8)}" y="${legendY - 5}" width="120" height="10" rx="2" `
        + `fill="url(#similarity-ramp)" stroke="${colors.border}" stroke-width="0.75"/>`);
    parts.push(`<text x="${r2(PADDING + dissimilarWidth + 136)}" y="${legendY}" dominant-baseline="central">`
        + 'Identical to the reference residue</text>');
    parts.push('</g>');

    parts.push('</svg>');
    return parts.join('\n');
}
