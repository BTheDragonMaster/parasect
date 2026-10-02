/**
 * A minimal ZIP writer, enough to bundle a handful of exported files.
 *
 * Pulling in a zip library for this would cost ~30 KB gzipped in a bundle that
 * already ships everything to every visitor, and the format's happy path is
 * small: local header, data, central directory, end record. Entries are
 * deflated with the browser's own CompressionStream where it exists (every
 * current browser) and stored uncompressed where it does not, which is still a
 * valid archive.
 */

const CRC_TABLE = (() => {
    const table = new Uint32Array(256);
    for (let i = 0; i < 256; i++) {
        let c = i;
        for (let k = 0; k < 8; k++) c = c & 1 ? 0xEDB88320 ^ (c >>> 1) : c >>> 1;
        table[i] = c >>> 0;
    }
    return table;
})();

function crc32(bytes) {
    let c = 0xFFFFFFFF;
    for (let i = 0; i < bytes.length; i++) c = CRC_TABLE[(c ^ bytes[i]) & 0xFF] ^ (c >>> 8);
    return (c ^ 0xFFFFFFFF) >>> 0;
}

async function deflateRaw(bytes) {
    if (typeof CompressionStream === 'undefined') return null;
    try {
        const stream = new Blob([bytes]).stream().pipeThrough(new CompressionStream('deflate-raw'));
        return new Uint8Array(await new Response(stream).arrayBuffer());
    } catch (err) {
        return null;  // unsupported algorithm: fall back to storing
    }
}

/** MS-DOS date/time, which is what the format stores. */
function dosDateTime(date) {
    const time = (date.getHours() << 11) | (date.getMinutes() << 5) | (Math.floor(date.getSeconds() / 2));
    const day = ((date.getFullYear() - 1980) << 9) | ((date.getMonth() + 1) << 5) | date.getDate();
    return { time, day };
}

class ByteWriter {
    constructor() {
        this.parts = [];
        this.length = 0;
    }

    u16(v) { this.raw(new Uint8Array([v & 0xFF, (v >>> 8) & 0xFF])); }

    u32(v) { this.raw(new Uint8Array([v & 0xFF, (v >>> 8) & 0xFF, (v >>> 16) & 0xFF, (v >>> 24) & 0xFF])); }

    raw(bytes) { this.parts.push(bytes); this.length += bytes.length; }
}

/**
 * Build a ZIP archive.
 *
 * @param {Array<{name: string, data: (Uint8Array|string)}>} files - entries to store.
 * @returns {Promise<Blob>} - the archive, ready to hand to a download link.
 */
export async function createZip(files) {
    const encoder = new TextEncoder();
    const { time, day } = dosDateTime(new Date());

    const out = new ByteWriter();
    const central = [];

    for (const file of files) {
        const name = encoder.encode(file.name);
        const raw = typeof file.data === 'string' ? encoder.encode(file.data) : file.data;
        const crc = crc32(raw);

        const deflated = await deflateRaw(raw);
        const useDeflate = deflated !== null && deflated.length < raw.length;
        const body = useDeflate ? deflated : raw;
        const method = useDeflate ? 8 : 0;

        const offset = out.length;
        out.u32(0x04034B50);           // local file header
        out.u16(20);                   // version needed
        out.u16(0x0800);               // UTF-8 names
        out.u16(method);
        out.u16(time);
        out.u16(day);
        out.u32(crc);
        out.u32(body.length);
        out.u32(raw.length);
        out.u16(name.length);
        out.u16(0);                    // no extra field
        out.raw(name);
        out.raw(body);

        const entry = new ByteWriter();
        entry.u32(0x02014B50);         // central directory header
        entry.u16(20);                 // version made by
        entry.u16(20);                 // version needed
        entry.u16(0x0800);
        entry.u16(method);
        entry.u16(time);
        entry.u16(day);
        entry.u32(crc);
        entry.u32(body.length);
        entry.u32(raw.length);
        entry.u16(name.length);
        entry.u16(0);                  // extra
        entry.u16(0);                  // comment
        entry.u16(0);                  // disk number
        entry.u16(0);                  // internal attrs
        entry.u32(0);                  // external attrs
        entry.u32(offset);
        entry.raw(name);
        central.push(entry);
    }

    const centralOffset = out.length;
    central.forEach((entry) => entry.parts.forEach((part) => out.raw(part)));
    const centralSize = out.length - centralOffset;

    out.u32(0x06054B50);               // end of central directory
    out.u16(0);                        // this disk
    out.u16(0);                        // disk with central directory
    out.u16(files.length);
    out.u16(files.length);
    out.u32(centralSize);
    out.u32(centralOffset);
    out.u16(0);                        // no archive comment

    return new Blob(out.parts, { type: 'application/zip' });
}

/**
 * Hand a blob to the browser as a download.
 *
 * @param {Blob} blob - file contents.
 * @param {string} filename - suggested name.
 * @returns {void}
 */
export function downloadBlob(blob, filename) {
    const url = URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = url;
    link.download = filename;
    document.body.appendChild(link);
    link.click();
    link.remove();
    // the object URL keeps the blob alive; release it once the click is handled
    setTimeout(() => URL.revokeObjectURL(url), 1000);
}
