from sys import argv
from paras.scripts.parsers.fasta import read_fasta


def sequences_are_equivalent(seq1: str, seq2: str, min_overlap_ratio: float = 0.8) -> bool:
    """
    Check if two DNA sequences represent the same sequence with partial overlaps.

    Parameters:
    - seq1, seq2: DNA sequences (strings)
    - min_overlap_ratio: minimum fraction of the shorter sequence length to consider as valid overlap

    Returns:
    - True if sequences overlap sufficiently and match in the overlapping region.
    - False otherwise.
    """

    def overlap_length(a, b):
        max_overlap = min(len(a), len(b))
        for length in range(max_overlap, 0, -1):
            if a[-length:] == b[:length]:
                return length
        return 0

    overlap_1_2 = overlap_length(seq1, seq2)
    overlap_2_1 = overlap_length(seq2, seq1)

    min_len = min(len(seq1), len(seq2))
    threshold = int(min_len * min_overlap_ratio)

    if overlap_1_2 >= threshold:
        return True
    if overlap_2_1 >= threshold:
        return True
    if seq1 in seq2:
        return True
    if seq2 in seq1:
        return True

    return False


def match_domains(domains_old, domains_new):
    id_to_seq_old = read_fasta(domains_old)
    id_to_seq_new = read_fasta(domains_new)

    old_matching = []
    old_mismatching = []
    old_no_match = []
    matching_ids = []
    mismatching_ids = []
    missing_ids = []

    for seq_id_1, seq_1 in id_to_seq_old.items():
        seq_ids = seq_id_1.split('|')
        match_found = False
        for seq_id in seq_ids:
            seq_id_found = False
            for seq_id_2, seq_2 in id_to_seq_new.items():

                if seq_id == seq_id_2:
                    seq_id_found = True
                    match_found = True
                    if sequences_are_equivalent(seq_1.upper(), seq_2.upper()):
                        old_matching.append(seq_id_1)
                        matching_ids.append(seq_id)

                    else:
                        old_mismatching.append(seq_id_1)
                        mismatching_ids.append(seq_id)
            if not seq_id_found:
                missing_ids.append(seq_id)

        if not match_found:
            old_no_match.append(seq_id_1)

    print(f"Matching: {len(old_matching)}")

    print(f"Mismatching: {len(old_mismatching)}")
    print(old_mismatching)
    partial_mismatches = list(set(old_matching).intersection(set(old_mismatching)))

    print(f"Partial mismatch: {len(partial_mismatches)}")
    print(partial_mismatches)

    print(f"Mismatching ids:")
    print(mismatching_ids)

    full_mismatches = list(set(old_mismatching) - set(old_matching))

    print(f"Full mismatches: {len(full_mismatches)}")
    print(full_mismatches)

    for partial_mismatch in partial_mismatches:
        ids = partial_mismatch.split('|')
        for seq_id in ids:
            if seq_id in mismatching_ids:
                print(seq_id)

    print(f"No match: {len(old_no_match)}")
    print(old_no_match)

    print(f"Missing ids: {len(missing_ids)}")
    print(missing_ids)

    return missing_ids


if __name__ == "__main__":
    match_domains(argv[1], argv[2])