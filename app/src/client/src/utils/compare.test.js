import {
    compareUrl, hammingDistance, moveItem, sortByDistance, parseCompareParams, residueSimilarity, similarityColors,
} from './compare';

const RAMP = { low: '#FFFFFF', high: '#1B5FA8', inks: ['#1C1D1A', '#FFFFFF'] };

describe('similarityColors', () => {
    it('runs from the low colour to the high colour', () => {
        expect(similarityColors(0, RAMP).fill).toBe('#ffffff');
        expect(similarityColors(1, RAMP).fill).toBe('#1b5fa8');
    });

    it('flips the ink as the fill darkens', () => {
        expect(similarityColors(0, RAMP).text).toBe('#1C1D1A');
        expect(similarityColors(1, RAMP).text).toBe('#FFFFFF');
    });

    it('clamps out-of-range values', () => {
        expect(similarityColors(-1, RAMP).fill).toBe('#ffffff');
        expect(similarityColors(2, RAMP).fill).toBe('#1b5fa8');
    });
});

describe('residueSimilarity', () => {
    const table = { alphabet: 'AL', matrix: [[1, 0.4], [0.4, 1]] };

    it('looks residues up case-insensitively', () => {
        expect(residueSimilarity(table, 'a', 'L')).toBe(0.4);
    });

    it('is null for unknown residues', () => {
        expect(residueSimilarity(table, 'A', 'Z')).toBeNull();
        expect(residueSimilarity(null, 'A', 'A')).toBeNull();
    });
});

describe('moveItem', () => {
    it('moves an item down and up', () => {
        expect(moveItem(['a', 'b', 'c'], 0, 2)).toEqual(['b', 'c', 'a']);
        expect(moveItem(['a', 'b', 'c'], 2, 0)).toEqual(['c', 'a', 'b']);
    });

    it('leaves the list alone for a no-op or bad index', () => {
        const list = ['a', 'b'];
        expect(moveItem(list, 1, 1)).toBe(list);
        expect(moveItem(list, 5, 0)).toBe(list);
    });
});

describe('compare URL', () => {
    it('round-trips ids and custom signatures, including a colon in a name', () => {
        const url = compareUrl({ ids: [3, 1, 2], custom: [{ name: 'my: query', signature: 'ABC' }] });
        const parsed = parseCompareParams(new URLSearchParams(url.split('?')[1]));
        expect(parsed).toEqual({ ids: [3, 1, 2], custom: [{ name: 'my: query', signature: 'ABC' }] });
    });

    it('ignores junk ids', () => {
        expect(parseCompareParams(new URLSearchParams('ids=1,x,,4')).ids).toEqual([1, 4]);
    });
});

describe('hammingDistance', () => {
    it('counts mismatching positions', () => {
        expect(hammingDistance('ABCD', 'ABXD')).toBe(1);
        expect(hammingDistance('ABCD', 'ABCD')).toBe(0);
    });

    it('counts missing positions as mismatches', () => {
        expect(hammingDistance('ABCD', 'AB')).toBe(2);
    });

    it('is null when either signature is missing', () => {
        expect(hammingDistance('', 'AB')).toBeNull();
    });
});

describe('sortByDistance', () => {
    const rows = [
        { key: 'ref', s: 'AAAA' },
        { key: 'far', s: 'XXXA' },
        { key: 'none', s: '' },
        { key: 'near', s: 'AAAX' },
        { key: 'near2', s: 'XAAA' },
    ];
    const keys = (list) => list.map((r) => r.key);

    it('keeps the reference on top and sorts the rest, missing last', () => {
        expect(keys(sortByDistance(rows, 's', 'asc'))).toEqual(['ref', 'near', 'near2', 'far', 'none']);
        expect(keys(sortByDistance(rows, 's', 'desc'))).toEqual(['ref', 'far', 'near', 'near2', 'none']);
    });
});
