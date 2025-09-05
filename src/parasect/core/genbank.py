# -*- coding: utf-8 -*-

"""Module for handling HMMER requests."""


import subprocess
from typing import List
import os
import itertools
from Bio import SeqIO

from ncbi_acc_download.errors import DownloadError

from parasect.core.writers import write_fasta_file


def genbank_to_fasta(path_in: str, path_out: str) -> None:
    """Parse protein sequences from a GenBank file and writes them to a fasta file.

    :param path_in: Path to input GenBank file.
    :type path_in: str
    :param path_out: Path to output fasta file.
    :type path_out: str
    :raises FileNotFoundError: If the file at the specified path does not exist.
    """
    # check if the file exists
    if not os.path.exists(path_in):
        raise FileNotFoundError(f"File not found: {path_in}")

    # initialize a counter for generating gene IDs
    counter = itertools.count()

    # initialize a dictionary to store the fasta sequences
    fasta_dict = {}

    # parse the GenBank file
    for record in SeqIO.parse(path_in, "genbank"):
        for feature in record.features:

            # check if the feature is a coding sequence (CDS)
            if feature.type == "CDS":

                # check if the feature has a translation
                if "translation" in feature.qualifiers:
                    sequence = feature.qualifiers["translation"]

                    # check if the feature has a protein ID, gene ID, or locus tag
                    # to use as the sequence ID. If not, generate a gene ID
                    if "protein_id" in feature.qualifiers:
                        seq_id = feature.qualifiers["protein_id"][0]
                    elif "gene_id" in feature.qualifiers:
                        seq_id = feature.qualifiers["gene_id"][0]
                    elif "locus_tag" in feature.qualifiers:
                        seq_id = feature.qualifiers["locus_tag"][0]
                    else:
                        seq_id = f"gene_{next(counter)}"

                    # add the sequence to the dictionary
                    fasta_dict[seq_id] = sequence

    # write the fasta sequences to a file
    write_fasta_file(fasta_dict=fasta_dict, path_out=path_out)


def fetch_from_genbank(protein_accessions: List[str], fasta_out: str) -> None:
    """Run hmmpfam2 on the given HMM database and fasta file.

    :param protein_accession: NCBI accession of protein sequence.
    :type protein_accession: str
    :param fasta_out: path to output fasta file
    :type fasta_out: str
    """

    command = ["ncbi-acc-download", "--format", "fasta", '--molecule', 'protein', '--out', fasta_out] + protein_accessions
    try:
        subprocess.check_call(command)
    except subprocess.CalledProcessError:
        raise DownloadError("Could not find one or more NCBI accessions")
