from paras.scripts.parsers.fasta import read_fasta, write_fasta
from argparse import ArgumentParser, Namespace
from paras.scripts.parsers.tabular import Tabular
import re
from typing import Optional, Any

from math import ceil
from collections import defaultdict
import os


def parse_args() -> Namespace:
    parser = ArgumentParser(description="Populate empty SQL database with data")
    parser.add_argument("--parasect", required=True, type=str, help="Path to parasect dataset file")
    parser.add_argument("--signature", required=True, type=str, help="Path to fasta containing A-domain signatures")
    parser.add_argument("--extended", required=True, type=str,
                        help="Path to fasta containing A-domain extended signatures")
    parser.add_argument("--domains", required=True, type=str,
                        help="Path to fasta containing A-domain sequences")
    parser.add_argument("--protein", required=True, type=str, help="Path to fasta file containing full protein sequences")
    parser.add_argument("--out", required=True, type=str, help="Path to output directory")
    args = parser.parse_args()

    return args

def write_tabular(dictionaries: list[dict[str, Any]], header: list[str], out_file: str) -> None:
    """
    Write tabular file from list of dictionaries, sorted by the first column

    :param dictionaries: list of dictionaries. Keys of all dictionaries should be the same
    :type dictionaries: list[dict[str, Any]]
    :param header: list of categories for labelling the header. Length must be one more than the list of dictionaries.
    First element represents category describing dictionary keys
    :type header: list[str]
    :param out_file: path to output file
    :type out_file: str
    """

    assert len(dictionaries) + 1 == len(header)
    assert len(header) >= 2
    assert len(dictionaries) >= 1

    with open(out_file, 'w') as out:

        header_string = '\t'.join(header)
        out.write(f"{header_string}\n")

        keys = list(dictionaries[0].keys())
        keys.sort()

        for key in keys:
            row = [key]
            for dictionary in dictionaries:
                row.append(dictionary[key])

            row_string = '\t'.join(row)
            out.write(f"{row_string}\n")


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


def deduplicate_sequences(id_to_seq: dict[str, str],
                          id_to_spec) -> tuple[dict[str, str], dict[str, str]]:
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
            if sequences_are_equivalent(id_to_seq[id1], id_to_seq[id2], min_overlap_ratio=0.8):

                union(id1, id2)
                if id_to_spec[id1] != id_to_spec[id2]:
                    print(f"MISMATCHING SPEC: {id1}, {id2}")

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


def remap_sequences(id_to_seq: dict[str, str], mapping: dict[str, str]) -> dict[str, str]:
    """
    Remap sequences to new identifiers

    :param id_to_seq: dictionary of old id to sequence
    :type id_to_seq: dict[str, str]
    :param mapping: dictionary of old id to new id
    :type mapping: dict[str, str]

    :return: dictionary of new id to sequence
    :rtype: dict[str, str]
    """

    new_id_to_seq: dict[str, str] = {}

    for old_id, seq in id_to_seq.items():
        new_id = mapping[old_id]
        if new_id in new_id_to_seq:
            if seq != new_id_to_seq[new_id]:
                raise ValueError(f"Mismatching signature sequences found for {new_id}: {seq}, {new_id_to_seq[new_id]}. Fix before updating database.")
        else:
            new_id_to_seq[new_id] = seq

    return new_id_to_seq


def deduplicate_domains(domain_mappings: list[dict[str, str]],
                        parasect_datasets: list[str]) -> tuple[dict[str, str], dict[str, str], dict[str, str]]:
    """Deduplicate domains by name and by sequence

    :param domain_mappings: list of dictionaries of sequence IDs to domain sequence, with each dictionary containing
                            the domains from a single user submission
    :type domain_mappings: list[dict[str, str]
    :param parasect_datasets: list of paths to parasect data files
    :type parasect_datasets: list[str]
    :return: dictionary of sequence ID to domain sequence, dictionary of old ids to new ids, dictionary of new ids to
    specificity
    :rtype: tuple[dict[str, str], dict[str, str], dict[str, str]]
    """

    id_to_seq: dict[str, str] = {}
    id_to_spec: dict[str, str] = {}

    annotated = 0

    for i, domain_mapping in enumerate(domain_mappings):
        parasect_data = Tabular(parasect_datasets[i], [0])
        for datapoint in parasect_data.data:
            seq_id = parasect_data.get_value(datapoint, 'domain_id')
            seq = domain_mapping[seq_id]

            specificity = parasect_data.get_value(datapoint, 'specificity')
            annotated += 1
            if seq_id in id_to_seq:
                if seq != id_to_seq[seq_id]:
                    raise ValueError(f"User-submitted domains of the same name found with mismatching sequences: {seq_id}. Please fix before updating database")
                if id_to_spec[seq_id] != specificity:
                    raise ValueError(
                        f"Inconsistent specificities reported for domain {seq_id}: {id_to_spec[seq_id]}, {specificity}. Fix before updating database.")
            else:
                id_to_seq[seq_id] = seq
                id_to_spec[seq_id] = specificity

    deduplicated_sequences, name_mapping = deduplicate_sequences(id_to_seq, id_to_spec)

    id_to_spec = remap_sequences(id_to_spec, name_mapping)

    print(f"Total domains: {annotated}")
    print(f"Total unique: {len(deduplicated_sequences)}")

    return deduplicated_sequences, name_mapping, id_to_spec


def deduplicate(domain_fasta, parasect_dataset, signature_file, extended_file, out_dir):
    domain_to_sig = read_fasta(signature_file)
    domain_to_ext = read_fasta(extended_file)
    domain_to_seq = read_fasta(domain_fasta)

    deduplicated_domains, mapping, id_to_spec = deduplicate_domains([domain_to_seq], [parasect_dataset])

    new_domain_to_sig = remap_sequences(domain_to_sig, mapping)
    new_domain_to_ext = remap_sequences(domain_to_ext, mapping)

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    new_domains = os.path.join(out_dir, "domains.fasta")
    write_fasta(deduplicated_domains, new_domains)
    new_sigs = os.path.join(out_dir, "signatures.fasta")
    write_fasta(new_domain_to_sig, new_sigs)
    new_ext = os.path.join(out_dir, "extended_signatures.fasta")
    write_fasta(new_domain_to_ext, new_ext)

    new_data = os.path.join(out_dir, "parasect_dataset.txt")
    write_tabular([deduplicated_domains, id_to_spec], ["domain_id", "sequence", "specificity"], new_data)


def map_domains_to_proteins(domain_fasta, protein_fasta, parasect_dataset, out_path):
    parasect_data = Tabular(parasect_dataset, [0])
    protein_to_seq = read_fasta(protein_fasta)
    counter = 0

    domain_to_seq = read_fasta(domain_fasta)
    domain_ids = list(domain_to_seq.keys())
    domain_ids.sort()

    for i, domain_1 in enumerate(domain_ids):
        seq_1 = domain_to_seq[domain_1]
        for domain_2 in domain_ids[i + 1:]:
            seq_2 = domain_to_seq[domain_2]
            if seq_1 == seq_2:
                print(f"Duplicate seqs found for {domain_1} and {domain_2}")

    used_proteins = set()

    for parasect_datapoint in parasect_data.data:
        domain = parasect_data.get_value(parasect_datapoint, "domain_id")
        seq = parasect_data.get_value(parasect_datapoint, "sequence")
        if seq != domain_to_seq[domain]:
            print(f"Mismatching sequences found for {domain}")

        domain_synonyms = domain.split('|')
        protein_found = False
        for domain_synonym in domain_synonyms:
            protein_synonym = '.'.join(domain_synonym.split('.')[:-1])
            if protein_synonym in protein_to_seq:
                protein_seq = protein_to_seq[protein_synonym]
                start_coords = [match.start() for match in re.finditer(seq, protein_seq)]

                end_coords = [start + len(seq) for start in start_coords]

                for start_coord in start_coords:
                    start_coord += 1

                if start_coords:
                    protein_found = True
                    used_proteins.add(protein_synonym)

                if len(start_coords) > 1:
                    print(f"Multiple matches found for {domain} in {protein_synonym}: {start_coords}")
                elif not start_coords:
                    print(f"No matches found for {domain} in {protein_synonym}")
                # else:
                #     print(start_coords[0], end_coords[0])

        if not protein_found:
            print(f"No protein found for {domain}")

    print(f"Found {counter} mismatching sequences")

    used_proteins = list(used_proteins)
    used_proteins.sort()

    with open(out_path, 'w') as out:
        for protein in used_proteins:
            out.write(f">{protein}\n{protein_to_seq[protein]}\n")



if __name__ == "__main__":
    args = parse_args()
    map_domains_to_proteins(args.domains, args.protein, args.parasect, args.out)