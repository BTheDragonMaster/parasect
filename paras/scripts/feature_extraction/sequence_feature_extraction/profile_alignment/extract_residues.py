#!/usr/bin/env python

from typing import List

import os
from sys import argv

from paras.scripts.parsers.fasta import read_fasta, write_fasta
from paras.scripts.feature_extraction.sequence_feature_extraction.read_positions import \
    POSITIONS_SIGNATURE, POSITIONS_EXTENDED_SIGNATURE, get_reference_positions
from paras.scripts.feature_extraction.sequence_feature_extraction.profile_alignment.align \
    import align_adomain, ALIGNMENT_FILE, REF_SEQUENCE
from paras.scripts.feature_extraction.sequence_feature_extraction.sequence_labels import parse_domain_id
from paras.scripts.feature_extraction.sequence_feature_extraction.adenylation_domain import AdenylationDomain

import multiprocessing as mp

CPUs = mp.cpu_count()


def stach_to_challis(id_to_stach):
    id_to_challis = {}
    for seq_id, stach in id_to_stach.items():
        challis = stach[:8]
        id_to_challis[seq_id] = challis

    return id_to_challis


def extract_stach_code(sequence):
    seq_id = "DOMAIN_TO_QUERY"

    aligned_domain, aligned_reference = align_adomain(seq_id, sequence, ALIGNMENT_FILE)

    aligned_positions_stach = get_reference_positions(POSITIONS_SIGNATURE, aligned_reference)
    aligned_positions_34 = get_reference_positions(POSITIONS_EXTENDED_SIGNATURE, aligned_reference)

    stach = []
    for position in aligned_positions_stach:
        stach.append(aligned_domain[position])
    stach = ''.join(stach)

    active_site = []
    for position in aligned_positions_34:
        active_site.append(aligned_domain[position])
    active_site = ''.join(active_site)
    return stach, active_site


def extract_stach_codes(fasta, alignment_file):

    id_to_seq = read_fasta(fasta)

    id_to_stach = {}
    id_to_34 = {}

    for seq_id, seq in id_to_seq.items():
        aligned_domain, aligned_reference = align_adomain(seq_id, seq, alignment_file)

        aligned_positions_stach = get_reference_positions(POSITIONS_SIGNATURE, aligned_reference)
        aligned_positions_34 = get_reference_positions(POSITIONS_EXTENDED_SIGNATURE, aligned_reference)

        stach = []

        for position in aligned_positions_stach:
            stach.append(aligned_domain[position])

        stach = ''.join(stach)
        id_to_stach[seq_id] = stach

        active_site = []
        for position in aligned_positions_34:
            active_site.append(aligned_domain[position])
        active_site = ''.join(active_site)
        id_to_34[seq_id] = active_site

        print(f"Processed {seq_id}")

    return id_to_stach, id_to_34


def get_stach_aa_signature(reference_alignment: str, domain_alignment: str) -> str:
    """ Extract stachelhaus residues from A domains """

    # Count residues in ref sequence and put positions in list
    poslist = get_reference_positions(POSITIONS_SIGNATURE, reference_alignment)
    # Extract positions from query sequence
    query_sig_seq = extract(domain_alignment, poslist)
    # Add fixed lysine 517
    query_sig_seq += "K"

    return query_sig_seq


def extract(sequence: str, positions: List[int]) -> str:
    """ Extracts a signature from an aligned sequence based on the provided
        positions. Accounts for gaps by looking behind or, if behind is already
        in the position list, ahead.

        Arguments:
            sequence: the aligned sequence to extract a signature from
            positions: the list of positions within the sequence to use

        Returns:
            the extracted signature as a string
    """
    seq = []
    for position in positions:
        aa = sequence[position]
        if aa == "-":
            if position - 1 not in positions:
                aa = sequence[position - 1]
            elif position + 1 not in positions:
                aa = sequence[position + 1]
        seq.append(aa)
    return "".join(seq)


def write_tabular(id_to_signature, out_file):
    with open(out_file, 'w') as out:
        for id, signature in id_to_signature.items():
            out.write(f'{id}\t{signature}\n')


def domains_from_fasta(fasta_in, job_name='paras_run'):
    """
    Extract adomains from fasta file and write them to out_dir

    Input:
    fasta_in: str, file location of .fasta file containing aa sequences
    job_name: str, name of job

    Output:
    fasta_out: str, file location of .fasta file containing detected adomains

    """

    id_to_seq = read_fasta(fasta_in)

    for seq_id, seq in id_to_seq.items():
        a_domain = AdenylationDomain()
        aligned_domain, aligned_reference = align_adomain(seq_id, seq, ALIGNMENT_FILE)

        aligned_positions_stach = get_reference_positions(POSITIONS_SIGNATURE, aligned_reference)
        aligned_positions_34 = get_reference_positions(POSITIONS_EXTENDED_SIGNATURE, aligned_reference)

        stach = []

        for position in aligned_positions_stach:
            stach.append(aligned_domain[position])

        stach = ''.join(stach)
        id_to_stach[seq_id] = stach

        active_site = []
        for position in aligned_positions_34:
            active_site.append(aligned_domain[position])
        active_site = ''.join(active_site)
        id_to_34[seq_id] = active_site

        print(f"Processed {seq_id}")

    return id_to_stach, id_to_34


if __name__ == "__main__":
    fasta = argv[1]
    alignment = ALIGNMENT_FILE
    out = argv[2]

    id_to_stach, id_to_34 = extract_stach_codes(fasta, alignment)

    if not os.path.exists(out):
        os.mkdir(out)

    write_fasta(id_to_stach, os.path.join(out, 'stachelhaus.fasta'))
    write_fasta(id_to_34, os.path.join(out, 'active_site_34.fasta'))
