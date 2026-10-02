import { checkSignature, describePlacement, parseSignatureInput } from './signatures';
import { makeDelimited } from './tabular';
import { NEIGHBOR_COLUMNS, neighborRows, queryRows } from './neighborExport';

const SIG_A = 'LDASFDASLFEVWGALTAGACLVLPPEEARKDPE';
const SIG_B = 'GWRSFP--WTKTKVFGGGEHNSYGPAEISIASHK';

describe('checkSignature', () => {
    it('accepts 34 residues regardless of case and whitespace', () => {
        expect(checkSignature(` ${SIG_A.slice(0, 10).toLowerCase()} ${SIG_A.slice(10)}\n`))
            .toEqual({ signature: SIG_A, problem: null });
    });

    it('accepts gaps and X', () => {
        expect(checkSignature(SIG_B).problem).toBeNull();
        expect(checkSignature(`X${SIG_A.slice(1)}`).problem).toBeNull();
    });

    it('names the characters that are not residues', () => {
        expect(checkSignature(`B${SIG_A.slice(1, 33)}Z`).problem).toBe("'B', 'Z' are not amino acids");
    });

    it('points out the 10-residue Stachelhaus code', () => {
        expect(checkSignature('DAWTIAAVCK').problem).toMatch(/10 residues, expected 34\. That looks like the 10-residue/);
    });

    it('points out a whole sequence', () => {
        expect(checkSignature(SIG_A.repeat(3)).problem).toMatch(/102 residues.*whole sequence/);
    });

    it('rejects empty input', () => {
        expect(checkSignature('   ').problem).toBe('no residues');
    });
});

describe('parseSignatureInput', () => {
    it('reads FASTA, joining wrapped sequence lines', () => {
        const text = `>dptA_A1 first\n${SIG_A.slice(0, 20)}\n${SIG_A.slice(20)}\n\n>\n${SIG_B}\n`;
        expect(parseSignatureInput(text)).toEqual({
            queries: [{ name: 'dptA_A1 first', signature: SIG_A }, { name: null, signature: SIG_B }],
            errors: [],
        });
    });

    it('reads bare signatures one per line', () => {
        expect(parseSignatureInput(`${SIG_A}\n\n${SIG_B.toLowerCase()}`).queries)
            .toEqual([{ name: null, signature: SIG_A }, { name: null, signature: SIG_B }]);
    });

    it('reads a name in front of a signature, as copied from a spreadsheet', () => {
        expect(parseSignatureInput(`dptA A1\t${SIG_A}\nsecond,${SIG_B}`).queries)
            .toEqual([{ name: 'dptA A1', signature: SIG_A }, { name: 'second', signature: SIG_B }]);
    });

    it('splits a line of several signatures', () => {
        expect(parseSignatureInput(`${SIG_A} ${SIG_B}`).queries).toHaveLength(2);
    });

    it('joins a signature typed with spaces in it', () => {
        expect(parseSignatureInput(`${SIG_A.slice(0, 17)} ${SIG_A.slice(17)}`).queries)
            .toEqual([{ name: null, signature: SIG_A }]);
    });

    it('reports each bad record by where it is', () => {
        const { queries, errors } = parseSignatureInput(`>good\n${SIG_A}\n>short\nDAWTIAAVCK`);
        expect(queries).toEqual([{ name: 'good', signature: SIG_A }]);
        expect(errors).toHaveLength(1);
        expect(errors[0]).toMatch(/^>short: 10 residues/);

        expect(parseSignatureInput(`${SIG_A}\nnot a signature`).errors[0]).toMatch(/^line 2: /);
    });
});

describe('describePlacement', () => {
    const cluster = (id, size) => ({ cluster_id: id, size, representative_name: `rep${id}` });

    it('says when a signature joins nothing', () => {
        expect(describePlacement({ clusters: [], within_threshold: 0, nearest_distance: 9 }, 5))
            .toMatch(/No reference domain is within 5 of 34 residues.*nearest is 9 away/);
    });

    it('names the one cluster a signature joins', () => {
        expect(describePlacement({ clusters: [cluster(3, 12)], within_threshold: 4, nearest_distance: 1 }, 5))
            .toBe('Joins the cluster around rep3 (12 domains): 4 reference domains within 5 residues, the nearest 1 away.');
    });

    it('says a signature bridging clusters would merge them', () => {
        const text = describePlacement(
            { clusters: [cluster(3, 12), cluster(8, 1)], within_threshold: 2, nearest_distance: 4 }, 5,
        );
        expect(text).toMatch(/^Bridges 2 clusters/);
        expect(text).toMatch(/rep8 \(a cluster of one\)/);
        expect(text).toMatch(/would merge them\.$/);
    });
});

describe('neighbour tables', () => {
    const neighbor = {
        rank: 1, id: 7, name: 'Q1.A1', distance: 3, substrates: ['2,3-diaminopropionic acid', 'serine'],
        genus: 'Streptomyces', kingdom: 'Bacteria', cluster_id: 42, extended_signature: SIG_B,
    };

    it('keeps multi-substrate cells in one TSV column', () => {
        const rows = neighborRows({ name: 'mine', signature: SIG_A, threshold: 5, neighbors: [neighbor] });
        const [header, line] = makeDelimited(NEIGHBOR_COLUMNS, rows, '\t').split('\n');
        expect(header.split('\t')).toHaveLength(NEIGHBOR_COLUMNS.length);
        const cells = line.split('\t');
        expect(cells).toHaveLength(NEIGHBOR_COLUMNS.length);
        expect(cells[header.split('\t').indexOf('substrates')]).toBe('2,3-diaminopropionic acid;serine');
        expect(cells[header.split('\t').indexOf('identity')]).toBe('0.912');
    });

    it('quotes CSV cells that contain the delimiter', () => {
        expect(makeDelimited([{ field: 'a' }], [{ a: 'x, "y"' }], ',')).toBe('a\n"x, ""y"""');
    });

    it('summarises each placed signature, skipping unplaced ones', () => {
        const queries = [{ key: 'q1', name: 'mine', signature: SIG_A }, { key: 'q2', name: 'pending', signature: SIG_B }];
        const placements = {
            q1: {
                nearest_distance: 3,
                within_threshold: 2,
                clusters: [{ cluster_id: 42, size: 10 }, { cluster_id: 7, size: 1 }],
                neighbors: [neighbor],
            },
        };
        expect(queryRows(queries, placements, 5)).toEqual([{
            query_name: 'mine',
            query_signature: SIG_A,
            threshold: 5,
            nearest_distance: 3,
            within_threshold: 2,
            clusters_joined: [42, 7],
            cluster_sizes: [10, 1],
            bridges_clusters: 'yes',
            nearest_neighbor: 'Q1.A1',
            nearest_neighbor_substrates: ['2,3-diaminopropionic acid', 'serine'],
        }]);
    });
});
