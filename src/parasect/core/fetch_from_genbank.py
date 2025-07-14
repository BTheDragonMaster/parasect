# -*- coding: utf-8 -*-

"""Module for handling HMMER requests."""


import subprocess
from typing import Dict, List, Tuple

from ncbi_acc_download.errors import DownloadError


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
