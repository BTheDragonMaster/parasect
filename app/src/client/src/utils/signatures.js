/**
 * Reading user-supplied extended signatures, and describing where they land in
 * the network. Mirrors the checks in server/routes/network.py (_clean_signature),
 * so a problem is reported per record before anything is sent.
 */

import { SIGNATURE_LENGTH } from '../components/DistanceBadge';

/** Same cap as the server's MAX_QUERIES. */
export const MAX_QUERIES = 50;
/** Same cap as the server's MAX_QUERY_NAME. */
export const MAX_QUERY_NAME = 120;

const RESIDUE = /[ACDEFGHIKLMNPQRSTVWYX-]/;
/** The classic Stachelhaus code, which is easy to paste by mistake. */
const SHORT_CODE_LENGTH = 10;

/**
 * Normalise one signature and say what, if anything, is wrong with it.
 *
 * @param {string} raw - the residues as pasted; whitespace and case don't matter.
 * @returns {{signature: string, problem: string|null}} - the cleaned signature, and a
 *     reason it can't be used, or null when it can.
 */
export function checkSignature(raw) {
    const signature = String(raw ?? '').replace(/\s+/g, '').toUpperCase();
    if (!signature) return { signature, problem: 'no residues' };

    const invalid = [...new Set([...signature].filter((c) => !RESIDUE.test(c)))];
    if (invalid.length) {
        const shown = invalid.slice(0, 5).map((c) => `'${c}'`).join(', ');
        return {
            signature,
            problem: `${shown} ${invalid.length === 1 ? 'is not an amino acid' : 'are not amino acids'}`,
        };
    }

    if (signature.length !== SIGNATURE_LENGTH) {
        let hint = '';
        if (signature.length === SHORT_CODE_LENGTH) {
            hint = '. That looks like the 10-residue Stachelhaus code, but this needs the 34-residue extended signature';
        } else if (signature.length > 2 * SIGNATURE_LENGTH) {
            hint = '. That looks like a whole sequence. Run it through PARAS first to get its extended signature';
        }
        return { signature, problem: `${signature.length} residues, expected ${SIGNATURE_LENGTH}${hint}` };
    }
    return { signature, problem: null };
}

function cleanName(name) {
    const trimmed = String(name ?? '').trim().slice(0, MAX_QUERY_NAME);
    return trimmed || null;
}

/**
 * Split pasted text into named records.
 *
 * FASTA if any line starts with '>'. Otherwise every non-empty line is one
 * record: either a signature on its own, or a name followed by a signature
 * (separated by whitespace, a comma or a semicolon), which is what a row copied
 * out of a spreadsheet looks like. A line of several complete signatures is
 * split into one record each.
 *
 * @param {string} text - the pasted input.
 * @returns {Array<{name: string|null, raw: string, where: string}>} - records in input order.
 */
function splitRecords(text) {
    const lines = String(text ?? '').split(/\r?\n/);
    const records = [];

    if (lines.some((line) => line.trim().startsWith('>'))) {
        let current = null;
        lines.forEach((line) => {
            const trimmed = line.trim();
            if (trimmed.startsWith('>')) {
                const name = cleanName(trimmed.slice(1));
                current = { name, raw: '', where: name ? `>${name}` : `record ${records.length + 1}` };
                records.push(current);
            } else if (trimmed) {
                if (!current) {
                    current = { name: null, raw: '', where: `record ${records.length + 1}` };
                    records.push(current);
                }
                current.raw += trimmed;
            }
        });
        return records;
    }

    lines.forEach((line, index) => {
        const trimmed = line.trim();
        if (!trimmed) return;
        const where = `line ${index + 1}`;
        const tokens = trimmed.split(/[\s,;]+/).filter(Boolean);
        const isSignature = (token) => !checkSignature(token).problem;

        if (tokens.length > 1 && tokens.every(isSignature)) {
            tokens.forEach((token) => records.push({ name: null, raw: token, where }));
        } else if (tokens.length > 1 && isSignature(tokens[tokens.length - 1])) {
            const last = tokens[tokens.length - 1];
            const name = cleanName(trimmed.slice(0, trimmed.lastIndexOf(last)).replace(/[\s,;]+$/, ''));
            records.push({ name, raw: last, where: name ? `"${name}"` : where });
        } else {
            // a signature typed with spaces in it, or something that isn't one:
            // either way the residue check below has the final say
            records.push({ name: null, raw: trimmed, where });
        }
    });
    return records;
}

/**
 * Parse pasted signatures.
 *
 * @param {string} text - FASTA, or one signature per line, optionally preceded by a name.
 * @returns {{queries: Array<{name: string|null, signature: string}>, errors: string[]}} -
 *     the usable records, and one message per record that isn't.
 */
export function parseSignatureInput(text) {
    const queries = [];
    const errors = [];
    splitRecords(text).forEach((record) => {
        const { signature, problem } = checkSignature(record.raw);
        if (problem) {
            errors.push(`${record.where}: ${problem}`);
        } else {
            queries.push({ name: record.name, signature });
        }
    });
    return { queries, errors };
}

/** "the cluster around Q04747.3.A2 (34 domains)". */
export function describeCluster(cluster) {
    const size = `${cluster.size} domain${cluster.size === 1 ? '' : 's'}`;
    return cluster.size === 1
        ? `${cluster.representative_name} (a cluster of one)`
        : `the cluster around ${cluster.representative_name} (${size})`;
}

/**
 * One sentence on where a placed signature landed.
 *
 * @param {object|undefined} placement - the server's placement for one query.
 * @param {number} threshold - the threshold it was placed at.
 * @returns {string} - the summary.
 */
export function describePlacement(placement, threshold) {
    if (!placement) return 'Placing...';
    const clusters = placement.clusters || [];
    const within = `${placement.within_threshold} reference domain${placement.within_threshold === 1 ? '' : 's'}`;

    if (!clusters.length) {
        return `No reference domain is within ${threshold} of 34 residues, so it joins no cluster. `
            + `The nearest is ${placement.nearest_distance} away, and its link is drawn fainter to show that.`;
    }
    if (clusters.length === 1) {
        return `Joins ${describeCluster(clusters[0])}: ${within} within ${threshold} residues, `
            + `the nearest ${placement.nearest_distance} away.`;
    }
    return `Bridges ${clusters.length} clusters at this threshold: ${within} within ${threshold} residues, `
        + `spread over ${clusters.map(describeCluster).join('; ')}. Added to the reference set, it would merge them.`;
}
