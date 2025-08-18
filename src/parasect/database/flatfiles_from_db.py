import os
from argparse import ArgumentParser, Namespace

from sqlalchemy import create_engine, select
from sqlalchemy.orm import Session

from parasect.database.build_database import AdenylationDomain, Protein, Substrate


def parse_arguments() -> Namespace:
    """ Parse command line arguments
    :return: arguments
    :rtype: Namespace
    """
    parser = ArgumentParser(description="Retrieve flatfiles from database")
    parser.add_argument('-db', required=True, type=str, help="Path to database")
    parser.add_argument('-o', required=True, type=str, help="Path to output directory")

    args = parser.parse_args()
    return args


def db_to_flatfiles(database_path: str, out_dir: str) -> None:
    """Write database entries to machine-readable flatfiles (.fasta, .tsv, .txt)

    :param database_path: path to sqlite database
    :type database_path: str
    :param out_dir: path to output directory
    :type out_dir: str

    """

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    engine = create_engine(f"sqlite:///{database_path}")
    domains_path = os.path.join(out_dir, "domains.fasta")
    proteins_path = os.path.join(out_dir, "proteins.fasta")
    substrates_path = os.path.join(out_dir, "smiles.tsv")
    parasect_path = os.path.join(out_dir, "parasect_dataset.txt")
    signatures_path = os.path.join(out_dir, "signatures.fasta")
    extended_signatures_path = os.path.join(out_dir, "extended_signatures.fasta")

    with Session(engine) as session:
        domains: list[AdenylationDomain] = list(session.scalars(select(AdenylationDomain)).all())
        proteins: list[Protein] = list(session.scalars(select(Protein)).all())
        substrates: list[Substrate] = list(session.scalars(select(Substrate)).all())

        domains.sort(key=lambda x: x.synonyms[0].synonym)
        proteins.sort(key=lambda x: x.synonyms[0].synonym)
        substrates.sort(key=lambda x: x.name)

        with open(domains_path, 'w') as domains_out, open(parasect_path, 'w') as parasect_out,\
                open(signatures_path, 'w') as signatures_out, open(extended_signatures_path, 'w') as ext_out:
            parasect_out.write("domain_id\tsequence\tspecificity\n")
            for domain in domains:
                domain_id = '|'.join([s.synonym for s in domain.synonyms])
                substrate_name = '|'.join([s.name for s in domain.substrates])
                domains_out.write(f">{domain_id}\n{domain.sequence}\n")
                parasect_out.write(f"{domain_id}\t{domain.sequence}\t{substrate_name}\n")
                signatures_out.write(f">{domain_id}\n{domain.signature}\n")
                ext_out.write(f">{domain_id}\n{domain.extended_signature}\n")
        with open(proteins_path, 'w') as proteins_out:
            for protein in proteins:
                protein_id = '|'.join([s.synonym for s in protein.synonyms])
                proteins_out.write(f">{protein_id}\n{protein.sequence}\n")
        with open(substrates_path, 'w') as substrates_out:
            substrates_out.write("substrate\tsmiles\n")
            for substrate in substrates:
                substrates_out.write(f"{substrate.name}\t{substrate.smiles}\n")


if __name__ == "__main__":
    args = parse_arguments()
    db_to_flatfiles(args.db, args.o)
