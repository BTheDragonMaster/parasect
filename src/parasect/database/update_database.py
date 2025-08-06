import os

from parasect.core.parsing import iterate_over_dir, write_fasta_file, SubstrateData
from parasect.core.chem import is_same_molecule
from parasect.database.query_database import get_substrates_from_smiles


def deduplicate_new_smiles(substrates: list[SubstrateData]) -> list[SubstrateData]:
    """

    :param substrates: list of substrates, each comprising substrate name and substrate SMILES
    :type substrates: list[SubstrateData]
    :return: deduplicated list of substrates
    :rtype: list[SubstrateData]
    """



def collate_user_submissions(submission_dir: str, out_dir: str) -> None:
    """
    :param submission_dir: Path to directory containing user submissions
    :type submission_dir: str
    :param out_dir: Path to output directory
    :type out_dir: str
    """
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    domains_out = os.path.join(out_dir, "domains.fasta")
    signatures_out = os.path.join(out_dir, "signatures.fasta")
    extended_out = os.path.join(out_dir, "extended_signatures.fasta")
    proteins_out = os.path.join(out_dir, "proteins.fasta")
    parasect_out = os.path.join(out_dir, "parasect_data.txt")
    smiles_out = os.path.join(out_dir, "smiles.tsv")

    for _, user_submission_path in iterate_over_dir(submission_dir, get_dirs=True):
        smiles_path = os.path.join(user_submission_path, "smiles.tsv")
        for submission_type, submission in iterate_over_dir(user_submission_path, get_dirs=True):
            domains_fasta = os.path.join(user_submission_path, "domains.fasta")
            signatures_fasta = os.path.join(user_submission_path, "signatures.fasta")
            extended_fasta = os.path.join(user_submission_path, "extended_signatures.fasta")
            parasect_data = os.path.join(user_submission_path, "parasect_data.txt")
            proteins_fasta = os.path.join(user_submission_path, "proteins.fasta")




