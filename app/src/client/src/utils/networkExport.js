/**
 * Export the network graph as a figure.
 *
 * The SVG is drawn from the graphology instance rather than scraped off sigma's
 * WebGL canvases, so it comes out as real vector geometry, and is editable in
 * Illustrator or Inkscape, which is what a figure usually needs. The export
 * reflects whatever is in the graph at that moment, expanded clusters included.
 * The PNG is rasterised from the same SVG so the two always agree.
 */

const PADDING = 48;
const LEGEND_SWATCH = 13;
const LEGEND_ROW = 21;
const LEGEND_PAD = 16;
const LEGEND_WIDTH = 260;
const LABEL_SIZE = 9;
/** Rough advance width for the label font at LABEL_SIZE; only used for collision boxes. */
const LABEL_CHAR_WIDTH = LABEL_SIZE * 0.52;
const FONT = "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif";

const escapeXml = (value) => String(value)
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;');

/** Round to 2dp; the coordinates carry far more precision than a figure needs. */
const r2 = (n) => Math.round(n * 100) / 100;

/**
 * Render the current graph to a standalone SVG document.
 *
 * @param {object} options - export options.
 * @param {import('graphology').default} options.graph - the live graph.
 * @param {(node: string, attrs: object) => string} options.colorOf - node colour, matching the on-screen reducer.
 * @param {boolean} options.showLabels - whether to draw node labels.
 * @param {Array<{label: string, color: string, count: number|undefined}>} options.legend - legend entries.
 * @param {{title: string, subtitle: string}} options.caption - figure heading.
 * @param {{background: string, text: string, textSecondary: string, edge: string, border: string}} options.colors - surface colours.
 * @param {number} options.size - width of the plot area in px.
 * @returns {string} - the SVG source.
 */
export function graphToSvg({ graph, colorOf, showLabels, legend, caption, colors, size = 1400 }) {
    let minX = Infinity;
    let maxX = -Infinity;
    let minY = Infinity;
    let maxY = -Infinity;
    graph.forEachNode((node, attrs) => {
        minX = Math.min(minX, attrs.x);
        maxX = Math.max(maxX, attrs.x);
        minY = Math.min(minY, attrs.y);
        maxY = Math.max(maxY, attrs.y);
    });
    if (!Number.isFinite(minX)) {
        minX = 0; maxX = 1; minY = 0; maxY = 1;
    }

    // square the extent so the layout keeps its aspect ratio
    const spanX = Math.max(maxX - minX, 1e-6);
    const spanY = Math.max(maxY - minY, 1e-6);
    const span = Math.max(spanX, spanY);
    const centreX = (minX + maxX) / 2;
    const centreY = (minY + maxY) / 2;

    const plot = size;
    const scale = (plot - 2 * PADDING) / span;
    // graph y grows upward, SVG y grows downward
    const px = (x) => r2(plot / 2 + (x - centreX) * scale);
    const py = (y) => r2(plot / 2 - (y - centreY) * scale);

    const headerHeight = caption?.title ? 64 : 0;
    const legendHeight = legend.length ? LEGEND_PAD * 2 + 24 + legend.length * LEGEND_ROW : 0;
    const legendX = plot + 1;
    const totalWidth = plot + (legend.length ? LEGEND_WIDTH : 0);
    const totalHeight = headerHeight + Math.max(plot, legendHeight);

    const parts = [];
    parts.push(`<?xml version="1.0" encoding="UTF-8"?>`);
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

    parts.push(`<g transform="translate(0 ${headerHeight})">`);

    // edges first so nodes sit on top of them
    parts.push(`<g stroke="${colors.edge}" stroke-width="0.8" stroke-opacity="0.75">`);
    graph.forEachEdge((edge, attrs, source, target, sourceAttrs, targetAttrs) => {
        parts.push(`<line x1="${px(sourceAttrs.x)}" y1="${py(sourceAttrs.y)}" `
            + `x2="${px(targetAttrs.x)}" y2="${py(targetAttrs.y)}"/>`);
    });
    parts.push('</g>');

    const placeable = [];
    parts.push('<g>');
    graph.forEachNode((node, attrs) => {
        const cx = px(attrs.x);
        const cy = py(attrs.y);
        const radius = r2(Math.max(attrs.size || 3, 1.5));
        parts.push(`<circle cx="${cx}" cy="${cy}" r="${radius}" fill="${colorOf(node, attrs)}"/>`);
        if (showLabels && attrs.label) placeable.push({ cx, cy, radius, text: attrs.label });
    });
    parts.push('</g>');

    if (placeable.length) {
        // On screen sigma drops labels that would collide; without the same pass the
        // export turns into a solid band of overlapping text. Biggest nodes win,
        // which is also what a figure wants to call out.
        placeable.sort((a, b) => b.radius - a.radius);
        const placed = [];
        const kept = [];
        for (const item of placeable) {
            const width = item.text.length * LABEL_CHAR_WIDTH;
            const box = {
                x1: item.cx + item.radius + 3,
                y1: item.cy - LABEL_SIZE * 0.6,
                x2: item.cx + item.radius + 3 + width,
                y2: item.cy + LABEL_SIZE * 0.6,
            };
            const collides = placed.some((b) => !(box.x2 < b.x1 || box.x1 > b.x2 || box.y2 < b.y1 || box.y1 > b.y2));
            if (collides) continue;
            placed.push(box);
            kept.push(`<text x="${r2(box.x1)}" y="${r2(item.cy + LABEL_SIZE * 0.35)}" `
                + `font-size="${LABEL_SIZE}" fill="${colors.text}">${escapeXml(item.text)}</text>`);
        }
        parts.push(`<g>${kept.join('')}</g>`);
    }

    if (legend.length) {
        parts.push(`<g transform="translate(${legendX} 0)">`);
        parts.push(`<line x1="0" y1="${LEGEND_PAD}" x2="0" y2="${Math.max(plot, legendHeight) - LEGEND_PAD}" `
            + `stroke="${colors.border}" stroke-width="1"/>`);
        parts.push(`<text x="${LEGEND_PAD}" y="${LEGEND_PAD + 16}" font-size="13" font-weight="700" `
            + `fill="${colors.text}">Legend</text>`);
        legend.forEach((item, i) => {
            const y = LEGEND_PAD + 24 + i * LEGEND_ROW;
            parts.push(`<rect x="${LEGEND_PAD}" y="${y + 4}" width="${LEGEND_SWATCH}" height="${LEGEND_SWATCH}" `
                + `rx="3" fill="${item.color}"/>`);
            const text = item.count === undefined ? item.label : `${item.label} (${item.count})`;
            parts.push(`<text x="${LEGEND_PAD + LEGEND_SWATCH + 8}" y="${y + 15}" font-size="11.5" `
                + `fill="${colors.text}">${escapeXml(text)}</text>`);
        });
        parts.push('</g>');
    }

    parts.push('</g></svg>');
    return parts.join('\n');
}

/**
 * Rasterise an SVG document.
 *
 * @param {string} svg - SVG source.
 * @param {number} scale - pixel density multiplier.
 * @returns {Promise<Blob>} - a PNG blob.
 */
export function svgToPng(svg, scale = 2) {
    return new Promise((resolve, reject) => {
        const widthMatch = svg.match(/width="(\d+(?:\.\d+)?)"/);
        const heightMatch = svg.match(/height="(\d+(?:\.\d+)?)"/);
        const width = widthMatch ? parseFloat(widthMatch[1]) : 1400;
        const height = heightMatch ? parseFloat(heightMatch[1]) : 1400;

        const image = new Image();
        // a data URL keeps the image same-origin, so the canvas stays untainted
        const url = `data:image/svg+xml;charset=utf-8,${encodeURIComponent(svg)}`;
        image.onload = () => {
            const canvas = document.createElement('canvas');
            canvas.width = Math.round(width * scale);
            canvas.height = Math.round(height * scale);
            const ctx = canvas.getContext('2d');
            ctx.drawImage(image, 0, 0, canvas.width, canvas.height);
            canvas.toBlob((blob) => (blob ? resolve(blob) : reject(new Error('canvas produced no image'))), 'image/png');
        };
        image.onerror = () => reject(new Error('could not rasterise the SVG'));
        image.src = url;
    });
}
