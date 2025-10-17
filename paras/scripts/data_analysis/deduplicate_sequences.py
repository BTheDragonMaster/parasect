from paras.scripts.parsers.fasta import read_fasta, write_fasta
from sys import argv
from math import ceil
from collections import defaultdict


def sequences_are_equivalent(seq1: str, seq2: str, min_overlap_ratio: float = 0.8) -> bool:
    seq1 = seq1.upper()
    seq2 = seq2.upper()

    def overlap_length(a, b):
        max_overlap = min(len(a), len(b))
        for length in range(max_overlap, 0, -1):
            if a[-length:] == b[:length]:
                return length
        return 0

    if seq1 == seq2:
        return True
    if seq1 in seq2 or seq2 in seq1:
        return True

    min_len = min(len(seq1), len(seq2))
    threshold = ceil(min_len * min_overlap_ratio)

    if overlap_length(seq1, seq2) >= threshold:
        return True
    if overlap_length(seq2, seq1) >= threshold:
        return True

    return False


def deduplicate_sequences(fasta_file, out_file, min_overlap_ratio: float = 0.8, choose_longest=True):
    id_to_seq = read_fasta(fasta_file)
    ids = sorted(id_to_seq.keys())

    # Union–find
    parent = {i: i for i in ids}
    rank = {i: 0 for i in ids}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra == rb:
            return
        if rank[ra] < rank[rb]:
            parent[ra] = rb
        elif rank[ra] > rank[rb]:
            parent[rb] = ra
        else:
            parent[rb] = ra
            rank[ra] += 1

    # Pairwise compare
    for i, id1 in enumerate(ids):
        for id2 in ids[i + 1:]:
            if sequences_are_equivalent(id_to_seq[id1], id_to_seq[id2], min_overlap_ratio):
                union(id1, id2)

    # Build groups
    groups = defaultdict(list)
    for i in ids:
        groups[find(i)].append(i)

    # Write output: join all IDs, pick representative
    dedup = {}
    for root, members in groups.items():
        members_sorted = sorted(members)
        if choose_longest:
            rep = max(members, key=lambda i: len(id_to_seq[i]))
        else:
            rep = members_sorted[0]
        new_id = "|".join(members_sorted)
        dedup[new_id] = id_to_seq[rep]

    write_fasta(dedup, out_file)


if __name__ == "__main__":
    deduplicate_sequences(argv[1], argv[2])