import os
from typing import TYPE_CHECKING

from paras.scripts.parsers.fasta import read_fasta

if TYPE_CHECKING:
    from paras.scripts.feature_extraction.sequence_feature_extraction.adenylation_domain import AdenylationDomain


MAPPING_SUFFIX: str = 'mapping.txt'
FASTA_SUFFIX: str = 'renamed_fasta.txt'


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

    mapping_file: str = os.path.join(out_dir, MAPPING_SUFFIX)
    new_fasta_file: str = os.path.join(out_dir, FASTA_SUFFIX)

    id_to_seq: dict[str, str] = read_fasta(fasta_file)
    counter: int = 0
    with open(new_fasta_file, 'w') as new_fasta:
        with open(mapping_file, 'w') as mapping:
            for seq_id, seq in id_to_seq.items():
                counter += 1
                mapping.write(f"{counter}\t{seq_id}\n")
                new_fasta.write(f">{counter}\n{seq}\n")

    return mapping_file, new_fasta_file


def parse_mapping(mapping_file):
    """
    Return a dictionary of renamed sequence IDs to original sequence IDs

    Parameters
    ----------
    mapping_file: str, path to mapping file which maps renamed sequence IDs to the original sequence IDs

    Returns
    -------
    new_to_original: dict of {new_id: old_id, ->}, with new_id and old_id strings

    """
    new_to_original = {}
    with open(mapping_file, 'r') as mapping:
        for line in mapping:
            line = line.strip()
            line_segments = line.split('\t')
            new = line_segments[0]
            original = '\t'.join(line_segments[1:])
            new_to_original[new] = original

    return new_to_original


def reverse_renaming(adenylation_domains: list["AdenylationDomain"], mapping_file: str) -> None:
    """
    Reverses the renaming of sequences within adenylation domain instances

    Parameters
    ----------

    adenylation_domains: list of [domain, ->], with each domain an AdenylationDomain instance
    mapping_file: str, path to mapping file which maps renamed sequence IDs to the original sequence IDs

    """
    new_to_original: dict[str, str] = parse_mapping(mapping_file)
    for domain in adenylation_domains:
        domain.protein_name = new_to_original[domain.protein_name]
