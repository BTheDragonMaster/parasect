# -*- coding: utf-8 -*-

"""Module for handling HMMER requests."""

import os
import subprocess
from typing import Dict, List, Tuple

from Bio import SearchIO
from Bio.SearchIO._model import HSP

from parasect.core.domain import AdenylationDomain
from parasect.core.parsing import parse_fasta_file


def run_hmmscan(hmm_dir, fasta_file, hmm_out):
    """
    Run hmmscan from command line

    Input:
    hmm_dir: str, dir of .hmm file containing the HMMs to be used in the scan
    fasta_dir: str, dir of .fasta file containing the sequences to be scanned
    out_dir: str, file location containing results of hmmscan

    """

    with open(hmm_out, 'w') as out:
        command = ['hmmscan', hmm_dir, fasta_file]
        subprocess.call(command, stdout=out)


def run_hmmpfam2(hmm_dir: str, fasta_file: str, hmm_out: str) -> None:
    """Run hmmpfam2 on the given HMM database and fasta file.

    :param hmm_dir: path to the HMM database file (.hmm).
    :type hmm_dir: str
    :param fasta_file: path to the fasta file.
    :type fasta_file: str
    :param hmm_out: path to the output file.
    :type hmm_out: str
    """
    with open(hmm_out, "w") as out:
        command = ["hmmpfam2", hmm_dir, fasta_file]
        subprocess.call(command, stdout=out)


def parse_hmm_results(path_in: str, hmmer_version: int = 2) -> Dict[str, HSP]:
    """Parse hmmpfam2 output file and return dictionary of domain identifier to Biopython HSP instance.

    :param path_in: path to hmmpfam2 output file (hmmer-2).
    :type path_in: str
    :param hmmer_version: version of HMMer
    "type hmmer_version: int, default: 2. Must be 2 or 3
    :return: Dictionary mapping domain identifier to Biopython HSP instance.
    :rtype: Dict[str, HSP]
    """

    if hmmer_version not in [2, 3]:
        raise ValueError(f"Unknown HMMer version: {hmmer_version}")
    filtered_hits = {}

    hmmer_string = f"hmmer{hmmer_version}-text"

    # parse relevant information from hmmpfam2 output
    for result in SearchIO.parse(path_in, hmmer_string):
        for hsp in result.hsps:

            # filter hits based on bitscore and hit_id
            if hsp.bitscore > 20:
                if hsp.hit_id == "AMP-binding" or hsp.hit_id == "AMP-binding_C":

                    header = f"{result.id}|{hsp.hit_id}|{hsp.query_start}-{hsp.query_end}"
                    filtered_hits[header] = hsp

    return filtered_hits


def rename_sequences(path_in: str, path_out: str) -> Tuple[str, str]:
    """Rename sequences before running hmmscan, and return file paths of mapping and new fasta file.

    :param path_in: path to input fasta file.
    :type path_in: str
    :param path_out: path to output directory.
    :type path_out: str
    :return: Tuple of mapping file path and new fasta file path. The mapping file
        maps renamed sequence IDs to the original sequence IDs.
    :rtype: Tuple[str, str]
    :raises FileNotFoundError: If the input file is not found.
    """
    # check if input file exists
    if not os.path.exists(path_in):
        raise FileNotFoundError(f"File {path_in} not found")

    # check if output exists, if not create it
    if not os.path.isdir(path_out):
        os.mkdir(path_out)

    # create paths for mapping and new fasta file
    path_mapping_file = os.path.join(path_out, "mapping.txt")
    path_new_fasta_file = os.path.join(path_out, "renamed_fasta.txt")

    # read fasta file
    id_to_seq = parse_fasta_file(path_in)

    counter = 0
    with open(path_new_fasta_file, "w") as fo_fasta:
        with open(path_mapping_file, "w") as fo_mapping:

            for seq_id, seq in id_to_seq.items():
                counter += 1

                fo_fasta.write(f">{counter}\n{seq}\n")
                fo_mapping.write(f"{counter}\t{seq_id}\n")

    return path_mapping_file, path_new_fasta_file


def reverse_renaming(
    adenylation_domains: List[AdenylationDomain], path_in_mapping_file: str
) -> None:
    """Reverse the renaming of sequences within adenylation domain instances.

    :param adenylation_domains: Adenylation domain instances.
    :type adenylation_domains: List[AdenylationDomain]
    :param path_in_mapping_file: path to mapping file which maps renamed sequence
        IDs to the original sequence IDs.
    :type path_in_mapping_file: str
    """
    new_to_original = {}

    with open(path_in_mapping_file, "r") as fo:
        for line in fo:
            line = line.strip()
            line_segments = line.split("\t")
            new = line_segments[0]
            original = "\t".join(line_segments[1:])
            new_to_original[new] = original

    for domain in adenylation_domains:
        domain.protein_name = new_to_original[domain.protein_name]
