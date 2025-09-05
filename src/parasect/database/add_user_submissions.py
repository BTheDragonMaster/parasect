import os
from argparse import ArgumentParser, Namespace
import traceback

from sqlalchemy import create_engine
from sqlalchemy.orm import Session

from parasect.core.parsing import iterate_over_dir, parse_fasta_file
from parasect.core.writers import write_fasta_file
from parasect.core.tabular import Tabular

from parasect.database.populate_database import populate_db

# TODO: Refactor such that entries are auto-added as pending to the database
# TODO: Create a table for proteins in the database
# TODO: Create a coordinate field for A-domains in the database


def parse_arguments() -> Namespace:
    """Parse arguments from command line

    :return: User-defined arguments
    :rtype: Namespace
    """
    parser = ArgumentParser(description="Add pending entries to database")
    parser.add_argument('-u', type=str, required=True, help="Path to user submissions folder")
    parser.add_argument('-db', type=str, required=True, help="Path to PARASECT database")

    args = parser.parse_args()
    return args


def replace_illegal_character(character: str, replacement_character: str, domains_file: str, signatures_file: str,
                              extended_file: str, proteins_file: str, dataset_file: str) -> None:
    """Replace illegal characters from sequence ids in files

    :param character: illegal character to replace
    :type character: str
    :param replacement_character: character to replace the illegal character with
    :type replacement_character: str
    :param domains_file: domains.fasta file path
    :type domains_file: str
    :param signatures_file: signatures.fasta file path
    :type signatures_file: str
    :param extended_file: extended_signatures.fasta file path
    :type extended_file: str
    :param proteins_file: proteins.fasta file path
    :type proteins_file: str
    :param dataset_file: parasect_data.txt file path
    :type dataset_file: str
    """
    protein_id_to_seq = parse_fasta_file(proteins_file)
    protein_name_mappings: dict[str, str] = {}
    new_protein_to_seq: dict[str, str] = {}

    for seq_id, protein_seq in protein_id_to_seq.items():
        new_seq_id = seq_id.replace(character, replacement_character)
        protein_name_mappings[seq_id] = new_seq_id
        new_protein_to_seq[new_seq_id] = protein_seq

    write_fasta_file(new_protein_to_seq, proteins_file)

    def get_domain_id(old_id: str) -> str:
        """ Get new domain ID from old domain ID
        :param old_id: old domain ID
        :type old_id: str
        :return: new domain ID
        :rtype: str
        """
        protein_id = '.'.join(old_id.split('.')[:-1])
        domain_number = old_id.split('.A')[-1]
        new_protein_id = protein_name_mappings[protein_id]
        new_domain_id = f"{new_protein_id}.A{domain_number}"
        return new_domain_id

    def change_domain_fasta_headers(fasta_file: str) -> None:
        """
        Change domain headers based on new protein names
        :param fasta_file: path to fasta file. File will be overwritten
        :type fasta_file: str
        """
        id_to_seq = parse_fasta_file(fasta_file)
        new_id_to_seq = {}
        for domain_id, seq in id_to_seq.items():
            new_id = get_domain_id(domain_id)
            new_id_to_seq[new_id] = seq

        write_fasta_file(new_id_to_seq, fasta_file)

    change_domain_fasta_headers(domains_file)
    change_domain_fasta_headers(signatures_file)
    change_domain_fasta_headers(extended_file)

    parasect_data = Tabular(dataset_file, separator='\t')

    with open(dataset_file, 'w') as parasect_out:
        parasect_out.write("domain_id\tsequence\tspecificity\n")
        for i, domain_name in enumerate(parasect_data.rows):
            new_domain_id = get_domain_id(domain_name)
            sequence = parasect_data.get_row_value(domain_name, "sequence")
            specificity = parasect_data.get_row_value(domain_name, "specificity")
            parasect_out.write(f"{new_domain_id}\t{sequence}\t{specificity}\n")


def remap_strings(id_to_string: dict[str, str], mapping: dict[str, str]) -> dict[str, str]:
    """
    Remap sequences to new identifiers

    :param id_to_string: dictionary of old id to sequence
    :type id_to_string: dict[str, str]
    :param mapping: dictionary of old id to new id
    :type mapping: dict[str, str]

    :return: dictionary of new id to sequence
    :rtype: dict[str, str]
    """

    new_id_to_string: dict[str, str] = {}

    for old_id, string in id_to_string.items():
        new_id = mapping[old_id]
        if new_id in new_id_to_string:
            if string != new_id_to_string[new_id]:
                raise ValueError(f"Mismatching strings found for {new_id}: {string}, {new_id_to_string[new_id]}. Fix before updating database.")
        else:
            new_id_to_string[new_id] = string

    return new_id_to_string

#TODO: Add correction handling

def add_user_submissions(submission_dir: str, database_path: str) -> None:
    """
    :param submission_dir: Path to directory containing user submissions
    :type submission_dir: str
    :param database_path: Path to parasect database
    :type database_path: str
    """
    engine = create_engine(f"sqlite:///{database_path}")

    with Session(engine) as session:

        try:
            for _, user_submission_path in iterate_over_dir(submission_dir, get_dirs=True):
                smiles_path = os.path.join(user_submission_path, "smiles.tsv")
                if not os.path.exists(smiles_path):
                    smiles_path = None
                for submission_type, submission in iterate_over_dir(user_submission_path, get_dirs=True):
                    if submission_type == 'new':
                        domains_fasta = os.path.join(submission, "domains.fasta")
                        signatures_fasta = os.path.join(submission, "signatures.fasta")
                        extended_fasta = os.path.join(submission, "extended_signatures.fasta")
                        parasect_data = os.path.join(submission, "parasect_data.txt")
                        proteins_fasta = os.path.join(submission, "proteins.fasta")

                        replace_illegal_character('|', '_', domains_fasta, signatures_fasta, extended_fasta,
                                                  proteins_fasta, parasect_data)

                        populate_db(session, parasect_data, smiles_path, signatures_fasta, extended_fasta,
                                    proteins_fasta)
                        session.flush()
            session.commit()
        except Exception as e:
            print(f"[ERROR] {type(e).__name__}: {e}")
            print("Rolling back changes.")
            traceback.print_exc()
            session.rollback()


def main() -> None:
    args = parse_arguments()
    add_user_submissions(args.u, args.db)


if __name__ == "__main__":
    main()
