from sys import argv
from collections import defaultdict

from paras.scripts.parsers.fasta import read_fasta, write_fasta

def deduplicate_sequences(id_to_seq: dict[str, str]) -> tuple[dict[str, str], dict[str, str]]:
    """
    Remove duplicate sequences, updating the sequence id

    :param id_to_seq: dictionary of sequence identifiers and sequences
    :type id_to_seq: dict[str, str]
    :param min_overlap_ratio: fraction of the shortest sequence (start and front) which has to overlap with the largest
    sequence in order for the sequences to be considered the same.
    :type min_overlap_ratio: float
    :return: dictionary of deduplicated sequences, with updated sequence IDs to reflect all original IDs; dictionary of
    old ids to new ids
    :rtype: tuple[dict[str, str], dict[str, str]]
    """

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
            if id_to_seq[id1] == id_to_seq[id2]:

                union(id1, id2)

    # Build groups
    groups = defaultdict(list)
    for i in ids:
        groups[find(i)].append(i)

    # Write output: join all IDs, pick representative
    dedup: dict[str, str] = {}
    old_to_new: dict[str, str] = {}
    for root, members in groups.items():
        members_sorted = sorted(members)
        rep = max(members, key=lambda i: len(id_to_seq[i]))

        new_id = "|".join(members_sorted)
        dedup[new_id] = id_to_seq[rep]
        for member in members_sorted:
            old_to_new[member] = new_id

    return dedup, old_to_new

def remove_duplicates(protein_file, protein_out):
    id_to_seq = read_fasta(protein_file)
    deduplicated, mapping = deduplicate_sequences(id_to_seq)
    protein_ids = list(deduplicated.keys())
    protein_ids.sort()
    with open(protein_out, 'w') as out:
        for protein_id in protein_ids:
            out.write(f">{protein_id}\n{deduplicated[protein_id]}\n")

if __name__ == "__main__":
    remove_duplicates(argv[1], argv[2])


    
