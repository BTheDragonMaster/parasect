/**
 * Draw an entity-relationship diagram of the database as a standalone SVG.
 */

const FONT = "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif";
const MONO = "ui-monospace, SFMono-Regular, Menlo, Consolas, monospace";

const PADDING = 32;
const CAPTION_HEIGHT = 58;
const LEGEND_HEIGHT = 34;
const HEADER_HEIGHT = 34;
const ROW_HEIGHT = 24;
const CELL_PAD = 12;
const BADGE_WIDTH = 24;
const COLUMN_GAP = 120;
const TABLE_GAP = 36;
const STUB = 18; // straight run out of a table, where the cardinality marks sit
const MIN_TABLE_WIDTH = 190;

const escapeXml = (value) => String(value)
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;');

const r2 = (n) => Math.round(n * 100) / 100;

let measureCtx = null;

/** Width of a string in a given CSS font; falls back to an estimate outside a browser. */
function textWidth(text, font) {
    if (measureCtx === null && typeof document !== 'undefined') {
        measureCtx = document.createElement('canvas').getContext('2d') || false;
    }
    if (!measureCtx) return String(text).length * 7.2;
    measureCtx.font = font;
    return measureCtx.measureText(String(text)).width;
}

const rowsLabel = (n) => `${n.toLocaleString('en-US')} ${n === 1 ? 'row' : 'rows'}`;

/** Tables' foreign keys as edges, ignoring self-references and dangling ones. */
function collectEdges(tables) {
    const byName = new Map(tables.map((t) => [t.name, t]));
    const edges = [];
    tables.forEach((table) => {
        table.columns.forEach((column, index) => {
            const ref = column.references;
            if (!ref || ref.table === table.name || !byName.has(ref.table)) return;
            const target = byName.get(ref.table);
            const targetIndex = Math.max(0, target.columns.findIndex((c) => c.name === ref.column));
            edges.push({
                from: table.name, fromRow: index, to: target.name, toRow: targetIndex, nullable: !column.notNull,
            });
        });
    });
    return edges;
}

/** Longest foreign-key path from each table to one without foreign keys. */
function computeLevels(tables, edges) {
    const out = new Map(tables.map((t) => [t.name, []]));
    edges.forEach((e) => out.get(e.from).push(e.to));

    const level = new Map();
    const visiting = new Set();
    const visit = (name) => {
        if (level.has(name)) return level.get(name);
        if (visiting.has(name)) return 0; // a cycle: break it here
        visiting.add(name);
        const value = out.get(name).reduce((max, next) => Math.max(max, visit(next) + 1), 0);
        visiting.delete(name);
        level.set(name, value);
        return value;
    };
    tables.forEach((t) => visit(t.name));
    return level;
}

/**
 * Render the schema to an SVG document.
 *
 * @param {object} options - render options.
 * @param {Array<{name: string, rowCount: number, columns: Array<{name: string, type: string, notNull: boolean,
 *     primaryKey: boolean, references: ({table: string, column: string}|null)}>}>} options.tables - schema.
 * @param {{title: string, subtitle: string}} options.caption - figure heading.
 * @param {{background: string, paper: string, header: string, border: string, text: string,
 *     textSecondary: string, edge: string, primaryKey: string, foreignKey: string}} options.colors - colours.
 * @returns {{svg: string, width: number, height: number}} - the SVG source and its size.
 */
export function schemaToSvg({ tables, caption, colors }) {
    const nameFont = `700 13px ${FONT}`;
    const metaFont = `400 11px ${FONT}`;
    const columnFont = `400 12px ${MONO}`;
    const typeFont = `400 11px ${MONO}`;

    const edges = collectEdges(tables);
    const levels = computeLevels(tables, edges);
    const maxLevel = Math.max(0, ...levels.values());

    // size every table from its own text
    const boxes = new Map(tables.map((table) => {
        const header = textWidth(table.name, nameFont) + 20 + textWidth(rowsLabel(table.rowCount), metaFont);
        const rows = table.columns.map((c) => BADGE_WIDTH + textWidth(c.name, columnFont) + 20 + textWidth(c.type, typeFont));
        const width = Math.ceil(Math.max(MIN_TABLE_WIDTH, header, ...rows) + 2 * CELL_PAD);
        const height = HEADER_HEIGHT + table.columns.length * ROW_HEIGHT;
        return [table.name, { table, width, height, x: 0, y: 0, column: maxLevel - levels.get(table.name) }];
    }));

    const columns = Array.from({ length: maxLevel + 1 }, () => []);
    [...boxes.values()]
        .sort((a, b) => a.table.name.localeCompare(b.table.name))
        .forEach((box) => columns[box.column].push(box));

    const columnWidths = columns.map((col) => Math.max(0, ...col.map((b) => b.width)));
    const columnX = [];
    columnWidths.reduce((x, w, i) => { columnX[i] = x; return x + w + COLUMN_GAP; }, PADDING);
    const columnHeight = (col) => col.reduce((h, b) => h + b.height, 0) + Math.max(0, col.length - 1) * TABLE_GAP;
    const contentHeight = Math.max(0, ...columns.map(columnHeight));
    const top = PADDING + CAPTION_HEIGHT;

    const place = () => {
        columns.forEach((col, i) => {
            let y = top + (contentHeight - columnHeight(col)) / 2;
            col.forEach((box) => {
                box.x = columnX[i];
                box.y = y;
                y += box.height + TABLE_GAP;
            });
        });
    };

    // order each column by where its neighbours sit, sweeping both ways a few times
    const neighbours = new Map([...boxes.keys()].map((name) => [name, []]));
    edges.forEach((e) => {
        neighbours.get(e.from).push(e.to);
        neighbours.get(e.to).push(e.from);
    });
    const centre = (box) => box.y + box.height / 2;
    place();
    for (let sweep = 0; sweep < 4; sweep++) {
        const order = sweep % 2 === 0 ? [...columns.keys()].reverse() : [...columns.keys()];
        order.forEach((i) => {
            const keyed = columns[i].map((box) => {
                const others = neighbours.get(box.table.name).map((n) => boxes.get(n)).filter((b) => b.column !== i);
                const key = others.length ? others.reduce((s, b) => s + centre(b), 0) / others.length : centre(box);
                return { box, key };
            });
            keyed.sort((a, b) => a.key - b.key);
            columns[i] = keyed.map((k) => k.box);
            place();
        });
    }

    const rowY = (box, row) => box.y + HEADER_HEIGHT + row * ROW_HEIGHT + ROW_HEIGHT / 2;

    /** Free vertical bands in a column: above, between and below its tables. */
    const gapsIn = (col) => {
        const gaps = [];
        let cursor = top - TABLE_GAP;
        col.forEach((box) => {
            gaps.push((cursor + box.y) / 2);
            cursor = box.y + box.height;
        });
        gaps.push(cursor + TABLE_GAP / 2);
        return gaps;
    };

    const edgePaths = edges.map((e) => {
        const src = boxes.get(e.from);
        const tgt = boxes.get(e.to);
        const sx = src.x + src.width;
        const sy = rowY(src, e.fromRow);
        const tx = tgt.x;
        const ty = rowY(tgt, e.toRow);

        // waypoints: out of the source, through a gap in every column in between, into the target
        const points = [[sx + STUB, sy]];
        for (let i = src.column + 1; i < tgt.column; i++) {
            const t = (columnX[i] - sx) / (tx - sx);
            const ideal = sy + (ty - sy) * t;
            const gy = gapsIn(columns[i]).reduce((best, g) => (Math.abs(g - ideal) < Math.abs(best - ideal) ? g : best));
            points.push([columnX[i] - 8, gy], [columnX[i] + columnWidths[i] + 8, gy]);
        }
        points.push([tx - STUB, ty]);

        let d = `M${r2(sx)},${r2(sy)} L${r2(points[0][0])},${r2(points[0][1])}`;
        for (let i = 1; i < points.length; i++) {
            const [x1, y1] = points[i - 1];
            const [x2, y2] = points[i];
            // odd steps cross a channel (curve), even steps cross a column (straight)
            if (i % 2 === 1) {
                const dx = (x2 - x1) / 2;
                d += ` C${r2(x1 + dx)},${r2(y1)} ${r2(x2 - dx)},${r2(y2)} ${r2(x2)},${r2(y2)}`;
            } else {
                d += ` L${r2(x2)},${r2(y2)}`;
            }
        }
        d += ` L${r2(tx)},${r2(ty)}`;

        // many: crow's foot opening onto the source table
        const many = `M${r2(sx + 12)},${r2(sy)} L${r2(sx)},${r2(sy - 6)} M${r2(sx + 12)},${r2(sy)} L${r2(sx)},${r2(sy + 6)}`;
        // one: double bar, or ring + bar for a nullable foreign key
        const one = e.nullable
            ? `M${r2(tx - 6)},${r2(ty - 6)} L${r2(tx - 6)},${r2(ty + 6)}`
            : `M${r2(tx - 6)},${r2(ty - 6)} L${r2(tx - 6)},${r2(ty + 6)} M${r2(tx - 11)},${r2(ty - 6)} L${r2(tx - 11)},${r2(ty + 6)}`;
        const ring = e.nullable
            ? `<circle cx="${r2(tx - 13)}" cy="${r2(ty)}" r="4" fill="${colors.background}" stroke="${colors.edge}" stroke-width="1.4"/>`
            : '';
        return `<g><path d="${d}" fill="none" stroke="${colors.edge}" stroke-width="1.4"/>`
            + `<path d="${many} ${one}" fill="none" stroke="${colors.edge}" stroke-width="1.4"/>${ring}</g>`;
    });

    const badge = (label, x, y, color) => `<text x="${r2(x)}" y="${r2(y)}" font-family="${FONT}" font-size="9"`
        + ` font-weight="700" fill="${color}" dominant-baseline="central">${label}</text>`;

    const tableGroups = [...boxes.values()].map((box) => {
        const { table, x, y, width, height } = box;
        const parts = [
            `<rect x="${r2(x)}" y="${r2(y)}" width="${width}" height="${height}" rx="8" fill="${colors.paper}" stroke="${colors.border}"/>`,
            // header: rounded on top only, so draw it clipped to the table outline
            `<path d="M${r2(x)},${r2(y + HEADER_HEIGHT)} V${r2(y + 8)} Q${r2(x)},${r2(y)} ${r2(x + 8)},${r2(y)}`
                + ` H${r2(x + width - 8)} Q${r2(x + width)},${r2(y)} ${r2(x + width)},${r2(y + 8)} V${r2(y + HEADER_HEIGHT)} Z"`
                + ` fill="${colors.header}"/>`,
            `<line x1="${r2(x)}" y1="${r2(y + HEADER_HEIGHT)}" x2="${r2(x + width)}" y2="${r2(y + HEADER_HEIGHT)}" stroke="${colors.border}"/>`,
            `<text x="${r2(x + CELL_PAD)}" y="${r2(y + HEADER_HEIGHT / 2)}" font-family="${FONT}" font-size="13" font-weight="700"`
                + ` fill="${colors.text}" dominant-baseline="central">${escapeXml(table.name)}</text>`,
            `<text x="${r2(x + width - CELL_PAD)}" y="${r2(y + HEADER_HEIGHT / 2)}" font-family="${FONT}" font-size="11"`
                + ` fill="${colors.textSecondary}" text-anchor="end" dominant-baseline="central">${rowsLabel(table.rowCount)}</text>`,
        ];
        table.columns.forEach((column, i) => {
            const cy = rowY(box, i);
            if (column.primaryKey) parts.push(badge('PK', x + CELL_PAD, cy - (column.references ? 5 : 0), colors.primaryKey));
            if (column.references) parts.push(badge('FK', x + CELL_PAD, cy + (column.primaryKey ? 5 : 0), colors.foreignKey));
            parts.push(
                `<text x="${r2(x + CELL_PAD + BADGE_WIDTH)}" y="${r2(cy)}" font-family="${MONO}" font-size="12"`
                    + ` font-weight="${column.primaryKey ? 700 : 400}" fill="${colors.text}" dominant-baseline="central">`
                    + `${escapeXml(column.name)}</text>`,
                `<text x="${r2(x + width - CELL_PAD)}" y="${r2(cy)}" font-family="${MONO}" font-size="11"`
                    + ` fill="${colors.textSecondary}" text-anchor="end" dominant-baseline="central">${escapeXml(column.type)}</text>`,
            );
        });
        return `<g>${parts.join('')}</g>`;
    });

    const lastColumn = columns.length - 1;
    const width = Math.ceil(columnX[lastColumn] + columnWidths[lastColumn] + PADDING);
    const legendY = top + contentHeight + 28;
    const height = Math.ceil(legendY + LEGEND_HEIGHT + PADDING - 20);

    // legend: key badges and the two cardinality marks
    const ly = legendY + 10;
    const legend = [
        badge('PK', PADDING, ly, colors.primaryKey),
        `<text x="${PADDING + 20}" y="${ly}" font-family="${FONT}" font-size="11" fill="${colors.textSecondary}" dominant-baseline="central">primary key</text>`,
        badge('FK', PADDING + 110, ly, colors.foreignKey),
        `<text x="${PADDING + 130}" y="${ly}" font-family="${FONT}" font-size="11" fill="${colors.textSecondary}" dominant-baseline="central">foreign key</text>`,
        `<path d="M${PADDING + 232},${ly} H${PADDING + 262} M${PADDING + 244},${ly} L${PADDING + 232},${ly - 6} M${PADDING + 244},${ly} L${PADDING + 232},${ly + 6}" fill="none" stroke="${colors.edge}" stroke-width="1.4"/>`,
        `<text x="${PADDING + 270}" y="${ly}" font-family="${FONT}" font-size="11" fill="${colors.textSecondary}" dominant-baseline="central">many</text>`,
        `<path d="M${PADDING + 322},${ly} H${PADDING + 352} M${PADDING + 341},${ly - 6} V${ly + 6} M${PADDING + 346},${ly - 6} V${ly + 6}" fill="none" stroke="${colors.edge}" stroke-width="1.4"/>`,
        `<text x="${PADDING + 360}" y="${ly}" font-family="${FONT}" font-size="11" fill="${colors.textSecondary}" dominant-baseline="central">exactly one</text>`,
    ];
    if (edges.some((e) => e.nullable)) {
        legend.push(
            `<path d="M${PADDING + 452},${ly} H${PADDING + 482} M${PADDING + 476},${ly - 6} V${ly + 6}" fill="none" stroke="${colors.edge}" stroke-width="1.4"/>`,
            `<circle cx="${PADDING + 469}" cy="${ly}" r="4" fill="${colors.background}" stroke="${colors.edge}" stroke-width="1.4"/>`,
            `<text x="${PADDING + 490}" y="${ly}" font-family="${FONT}" font-size="11" fill="${colors.textSecondary}" dominant-baseline="central">zero or one</text>`,
        );
    }

    const svg = [
        `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">`,
        `<rect width="${width}" height="${height}" fill="${colors.background}"/>`,
        `<text x="${PADDING}" y="${PADDING + 12}" font-family="${FONT}" font-size="18" font-weight="700" fill="${colors.text}">${escapeXml(caption.title)}</text>`,
        `<text x="${PADDING}" y="${PADDING + 34}" font-family="${FONT}" font-size="12" fill="${colors.textSecondary}">${escapeXml(caption.subtitle)}</text>`,
        ...edgePaths,
        ...tableGroups,
        ...legend,
        '</svg>',
    ].join('\n');

    return { svg, width, height };
}
