# -*- coding: utf-8 -*-

"""Module for handling HMMER requests."""

import os
import subprocess

from Bio import SearchIO
from Bio.SearchIO._model import HSP

from parasect.core.abstractions import AdenylationDomain
from parasect.core.parsing import read_fasta_file


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
        command = ["hmm2pfam", hmm_dir, fasta_file]
        subprocess.call(command, stdout=out)


def make_domain_id(seq_id: str, hit_id: str, start: int, end: int) -> str:
    """Compose a domain identifier from the sequence ID, hit ID, and start and end positions.

    :param seq_id: sequence ID.
    :type seq_id: str
    :param hit_id: hit ID.
    :type hit_id: str
    :param start: start position.
    :type start: int
    :param end: end position.
    :type end: int
    :return: domain identifier.
    :rtype: str

    .. note:: Uses '|' as a separator between the sequence ID, hit ID, and start 
        and end positions.
    """
    return f"{seq_id}|{hit_id}|{start}-{end}"


###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################


def parse_hmm2_results(hmm_results):
    """
    Return dictionary of domain identifier to Biopython HSP instance

    Parameters
    ----------
    hmm_results: str, path to hmmpfam2 output file (hmmer-2)

    Returns
    -------
    id_to_hit: dict of {domain_id: HSP, ->}, with domain_id str and HSP a Biopython HSP instance containing a HMM hit
        between an Hmm2 adenylation domain HMM and a query sequence

    """
    id_to_hit: dict[str, HSP] = {}

    for result in SearchIO.parse(hmm_results, "hmmer2-text"):
        for hsp in result.hsps:
            if hsp.bitscore > 20:
                if hsp.hit_id == "AMP-binding" or hsp.hit_id == "AMP-binding_C":
                    header = make_domain_id(result.id, hsp.hit_id, hsp.query_start, hsp.query_end)
                    id_to_hit[header] = hsp

    return id_to_hit


def rename_sequences(fasta_file: str, out_dir: str) -> tuple[str, str]:
    """

    Rename sequences before running hmmscan, and return file paths of mapping and new fasta file

    Parameters
    ----------
    fasta_file: str, path to input fasta file
    out_dir: str, path to output directory

    Returns
    -------
    mapping_file: str, path to mapping file which maps renamed sequence IDs to the original sequence IDs
    new_fasta_file: str, path to output fasta file

    """
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    mapping_file: str = os.path.join(out_dir, "mapping.txt")
    new_fasta_file: str = os.path.join(out_dir, "renamed_fasta.txt")

    id_to_seq: dict[str, str] = read_fasta_file(fasta_file)
    counter: int = 0
    with open(new_fasta_file, "w") as new_fasta:
        with open(mapping_file, "w") as mapping:
            for seq_id, seq in id_to_seq.items():
                counter += 1
                mapping.write(f"{counter}\t{seq_id}\n")
                new_fasta.write(f">{counter}\n{seq}\n")

    return mapping_file, new_fasta_file


def reverse_renaming(adenylation_domains: list["AdenylationDomain"], mapping_file: str) -> None:
    """
    Reverses the renaming of sequences within adenylation domain instances

    Parameters
    ----------

    adenylation_domains: list of [domain, ->], with each domain an AdenylationDomain instance
    mapping_file: str, path to mapping file which maps renamed sequence IDs to the original sequence IDs

    """
    new_to_original = {}
    with open(mapping_file, "r") as mapping:
        for line in mapping:
            line = line.strip()
            line_segments = line.split("\t")
            new = line_segments[0]
            original = "\t".join(line_segments[1:])
            new_to_original[new] = original

    for domain in adenylation_domains:
        domain.protein_name = new_to_original[domain.protein_name]
