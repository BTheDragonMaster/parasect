/**
 * Export the prediction results as a figure.
 */

import { FONT, escapeXml, fitText, r2, textWidth } from './svgText';

const PADDING = 32;
const AXIS_WIDTH = 40;
const HEADER_HEIGHT = 64;
const GENE_LABEL_HEIGHT = 24;
const PANEL_PAD = 10;
const VALUE_LABEL_HEIGHT = 16;
const PLOT_HEIGHT = 160;
const BAR_WIDTH = 22;
const BAR_GAP = 4;
const MIN_COLUMN_WIDTH = 56;
const DOMAIN_GAP = 14;
const GENE_GAP = 20;
const SUBSTRATE_LABEL_SIZE = 10;
/** Longer substrate names are cut short, so one odd name can't stretch the whole figure. */
const MAX_SUBSTRATE_LABEL_WIDTH = 140;
const DOMAIN_BOX_HEIGHT = 24;
const RANGE_LABEL_HEIGHT = 16;
const FOOTER_HEIGHT = 36;
const FOOTNOTE_HEIGHT = 18;

/**
 * Render the results to a standalone SVG document.
 *
 * @param {object} options - export options.
 * @param {Array<{proteinName: string, items: Array<object>}>} options.groups - protein groups in
 *     display order, each with its domain results in order along the protein.
 * @param {number} options.topN - how many predictions to plot per domain.
 * @param {{title: string, subtitle: string, footnote: (string|null|undefined)}} options.caption - figure
 *     heading, and an optional line of small print under the legend.
 * @param {{background: string, text: string, textSecondary: string, grid: string, topBar: string,
 *     otherBar: string, domainFill: string, domainText: string,
 *     bands: Array<{surface: string, rail: string}>}} options.colors - colours.
 * @returns {string} - the SVG source.
 */
export function resultsToSvg({ groups, topN, caption, colors }) {
    const columnWidth = Math.max(topN * BAR_WIDTH + (topN - 1) * BAR_GAP, MIN_COLUMN_WIDTH);

    // prepare every domain column once: its top predictions and their labels
    const layoutGroups = groups.map((group) => ({
        name: group.proteinName,
        domains: group.items.map((result) => ({
            result,
            // the top prediction's label is set in bold, so it is measured in bold too
            bars: (result.predictions || []).slice(0, topN).map((p, k) => ({
                name: p.substrate_name,
                probability: Math.min(Math.max(parseFloat(p.probability) || 0, 0), 1),
                weight: k === 0 ? 700 : 400,
                label: fitText(p.substrate_name, MAX_SUBSTRATE_LABEL_WIDTH, SUBSTRATE_LABEL_SIZE, k === 0 ? 700 : 400),
            })),
        })),
    }));

    const labelBand = Math.ceil(Math.max(
        20,
        ...layoutGroups.flatMap((g) => g.domains.flatMap((d) => d.bars.map(
            (b) => textWidth(b.label, SUBSTRATE_LABEL_SIZE, b.weight),
        ))),
    )) + 8;

    // vertical positions, relative to the top of a gene panel
    const plotTop = PANEL_PAD + VALUE_LABEL_HEIGHT;
    const baseline = plotTop + PLOT_HEIGHT;
    const domainTop = baseline + labelBand;
    const panelHeight = domainTop + DOMAIN_BOX_HEIGHT + RANGE_LABEL_HEIGHT + PANEL_PAD;

    const panelsTop = HEADER_HEIGHT + GENE_LABEL_HEIGHT;
    const y = (probability) => r2(panelsTop + baseline - probability * PLOT_HEIGHT);

    // horizontal positions
    let cursor = PADDING + AXIS_WIDTH;
    layoutGroups.forEach((group, i) => {
        if (i > 0) cursor += GENE_GAP;
        group.x = cursor;
        group.width = 2 * PANEL_PAD + group.domains.length * columnWidth
            + Math.max(group.domains.length - 1, 0) * DOMAIN_GAP;
        group.domains.forEach((domain, j) => {
            domain.x = cursor + PANEL_PAD + j * (columnWidth + DOMAIN_GAP);
        });
        cursor += group.width;
    });

    const plotLeft = PADDING + AXIS_WIDTH;
    const plotRight = cursor;
    const titleWidth = Math.max(
        textWidth(caption?.title || '', 19, 700),
        textWidth(caption?.subtitle || '', 12.5),
        textWidth(caption?.footnote || '', 10),
    );
    const totalWidth = Math.ceil(Math.max(plotRight + PADDING, PADDING * 2 + titleWidth));
    const totalHeight = panelsTop + panelHeight + FOOTER_HEIGHT + (caption?.footnote ? FOOTNOTE_HEIGHT : 0);

    const parts = [];
    parts.push('<?xml version="1.0" encoding="UTF-8"?>');
    parts.push(`<svg xmlns="http://www.w3.org/2000/svg" width="${totalWidth}" height="${totalHeight}" `
        + `viewBox="0 0 ${totalWidth} ${totalHeight}" font-family="${FONT}">`);
    parts.push(`<rect width="${totalWidth}" height="${totalHeight}" fill="${colors.background}"/>`);

    if (caption?.title) {
        parts.push(`<text x="${PADDING}" y="30" font-size="19" font-weight="700" fill="${colors.text}">`
            + `${escapeXml(caption.title)}</text>`);
        if (caption.subtitle) {
            parts.push(`<text x="${PADDING}" y="50" font-size="12.5" fill="${colors.textSecondary}">`
                + `${escapeXml(caption.subtitle)}</text>`);
        }
    }

    // gene panels: tinted band with a rail on the left and the gene's name above,
    // alternating the two band colours like on the results page
    layoutGroups.forEach((group, i) => {
        const band = colors.bands[i % colors.bands.length];
        const top = panelsTop;
        parts.push(`<g>`);
        parts.push(`<title>${escapeXml(group.name)}</title>`);
        parts.push(`<rect x="${group.x}" y="${top}" width="${group.width}" height="${panelHeight}" rx="10" `
            + `fill="${band.surface}"/>`);
        parts.push(`<rect x="${group.x}" y="${top}" width="4" height="${panelHeight}" fill="${band.rail}"/>`);
        parts.push(`<rect x="${group.x}" y="${top - 3}" width="${group.width}" height="3" fill="${band.rail}"/>`);
        const count = `${group.domains.length} domain${group.domains.length === 1 ? '' : 's'}`;
        const countWidth = textWidth(count, 10);
        // the name wins over the domain count: the count is dropped before the name is cut short
        const nameRoom = group.width - 4;
        const name = fitText(group.name, nameRoom, 12, 700);
        const nameWidth = textWidth(name, 12, 700);
        parts.push(`<text x="${group.x + 2}" y="${top - 9}" font-size="12" font-weight="700" `
            + `fill="${colors.text}">${escapeXml(name)}</text>`);
        if (name === group.name && nameWidth + 6 + countWidth <= nameRoom) {
            parts.push(`<text x="${r2(group.x + 2 + nameWidth + 6)}" y="${top - 9}" font-size="10" `
                + `fill="${colors.textSecondary}">${escapeXml(count)}</text>`);
        }
        parts.push('</g>');
    });

    // y axis and gridlines, drawn over the panels so they read across genes
    parts.push(`<g font-size="10" fill="${colors.textSecondary}">`);
    [0, 0.25, 0.5, 0.75, 1].forEach((tick) => {
        const ty = y(tick);
        parts.push(`<line x1="${plotLeft - 4}" y1="${ty}" x2="${plotRight}" y2="${ty}" stroke="${colors.grid}" `
            + `stroke-width="1"${tick === 0 ? '' : ' stroke-dasharray="3 3"'}/>`);
        parts.push(`<text x="${plotLeft - 8}" y="${ty}" text-anchor="end" dominant-baseline="central">`
            + `${tick.toFixed(2)}</text>`);
    });
    const axisMid = r2(panelsTop + plotTop + PLOT_HEIGHT / 2);
    parts.push(`<text transform="translate(${PADDING - 6} ${axisMid}) rotate(-90)" text-anchor="middle" `
        + `dominant-baseline="hanging" font-size="11">Probability</text>`);
    parts.push('</g>');

    // per domain: bars, value above, substrate name below, domain box
    layoutGroups.forEach((group) => {
        group.domains.forEach((domain) => {
            const { result, bars } = domain;
            const barsWidth = bars.length * BAR_WIDTH + Math.max(bars.length - 1, 0) * BAR_GAP;
            const barsLeft = domain.x + (columnWidth - barsWidth) / 2;
            parts.push('<g>');
            bars.forEach((bar, k) => {
                const bx = r2(barsLeft + k * (BAR_WIDTH + BAR_GAP));
                const cx = r2(bx + BAR_WIDTH / 2);
                const top = y(bar.probability);
                const height = r2(y(0) - top);
                const fill = k === 0 ? colors.topBar : colors.otherBar;
                parts.push(`<rect x="${bx}" y="${top}" width="${BAR_WIDTH}" height="${height}" rx="2" fill="${fill}">`
                    + `<title>${escapeXml(`${bar.name}: ${bar.probability.toFixed(3)}`)}</title></rect>`);
                parts.push(`<text x="${cx}" y="${r2(top - 4)}" text-anchor="middle" font-size="9" `
                    + `fill="${colors.textSecondary}">${bar.probability.toFixed(2)}</text>`);
                parts.push(`<text transform="translate(${cx} ${r2(y(0) + 6)}) rotate(-90)" text-anchor="end" `
                    + `dominant-baseline="central" font-size="${SUBSTRATE_LABEL_SIZE}" `
                    + `font-weight="${bar.weight}" fill="${colors.text}">${escapeXml(bar.label)}</text>`);
            });

            const boxY = panelsTop + domainTop;
            parts.push(`<rect x="${domain.x}" y="${boxY}" width="${columnWidth}" height="${DOMAIN_BOX_HEIGHT}" `
                + `rx="5" fill="${colors.domainFill}"/>`);
            parts.push(`<text x="${r2(domain.x + columnWidth / 2)}" y="${r2(boxY + DOMAIN_BOX_HEIGHT / 2)}" `
                + `text-anchor="middle" dominant-baseline="central" font-size="11.5" font-weight="700" `
                + `fill="${colors.domainText}">A${escapeXml(result.domain_nr)}</text>`);
            const range = `${result.domain_start}-${result.domain_end}`;
            parts.push(`<text x="${r2(domain.x + columnWidth / 2)}" y="${boxY + DOMAIN_BOX_HEIGHT + 12}" `
                + `text-anchor="middle" font-size="9" fill="${colors.textSecondary}">`
                + `${escapeXml(fitText(range, columnWidth, 9))}</text>`);
            parts.push('</g>');
        });
    });

    // footer: what the bars are
    const footerY = panelsTop + panelHeight + 22;
    parts.push(`<rect x="${PADDING}" y="${footerY - 9}" width="11" height="11" rx="2" fill="${colors.topBar}"/>`);
    parts.push(`<text x="${PADDING + 16}" y="${footerY}" font-size="11" fill="${colors.text}">Top prediction</text>`);
    const otherX = PADDING + 16 + textWidth('Top prediction', 11) + 18;
    if (topN > 1) {
        parts.push(`<rect x="${r2(otherX)}" y="${footerY - 9}" width="11" height="11" rx="2" fill="${colors.otherBar}"/>`);
        parts.push(`<text x="${r2(otherX + 16)}" y="${footerY}" font-size="11" fill="${colors.text}">`
            + `Next ${topN - 1} prediction${topN === 2 ? '' : 's'}</text>`);
    }

    if (caption?.footnote) {
        parts.push(`<text x="${PADDING}" y="${footerY + FOOTNOTE_HEIGHT + 2}" font-size="10" `
            + `fill="${colors.textSecondary}">${escapeXml(caption.footnote)}</text>`);
    }

    parts.push('</svg>');
    return parts.join('\n');
}
