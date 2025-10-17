from argparse import ArgumentParser, Namespace
import os

from sqlalchemy.orm import Session
from sqlalchemy import select, create_engine

from parasect.core.tabular import Tabular
from parasect.database.query_database import get_domains_from_protein_synonym
from parasect.database.build_database import AdenylationDomain
from parasect.core.writers import write_list


def parse_arguments() -> Namespace:
    """Parse arguments from command line

    :return: Arguments
    :rtype: Namespace
    """
    parser = ArgumentParser(description="Retrieve taxonomy for database domains")
    parser.add_argument("-db", '--database', required=True, type=str,
                        help="Path to PARASECT database")
    parser.add_argument("-t", '--taxonomy', required=True, type=str,
                        help="Path to taxonomy file, with clades listed in the header")
    parser.add_argument("-r", '--rank', type=str, default="family", help="Taxonomic rank to split on")
    parser.add_argument("-o", '--out', required=True, type=str,
                        help="Path to output directory")

    arguments = parser.parse_args()

    return arguments


def split_by_clade(session: Session, taxonomy_file: str, rank: str, out_dir: str) -> None:
    """Divide domains into files by clade

    :param session: database session
    :type session: Session
    :param taxonomy_file: path to file containing taxonomy, with taxonomic ranks indicated in the header
    :type taxonomy_file: str
    :param rank: taxonomic rank to split on
    :type rank: str
    :param out_dir: output directory
    :type out_dir: str
    """

    tax_data = Tabular(taxonomy_file)
    protein_to_rank: dict[str, str] = {}
    for protein in tax_data.rows:
        protein_to_rank[protein] = tax_data.get_row_value(protein, rank)

    rank_to_domains: dict[str, list[str]] = {}

    query = select(AdenylationDomain)
    all_domains = set(session.scalars(query))

    for protein, rank in protein_to_rank.items():
        domains = get_domains_from_protein_synonym(session, protein)
        if rank not in rank_to_domains:
            rank_to_domains[rank] = []
        for domain in domains:
            all_domains.discard(domain)
            rank_to_domains[rank].append(domain.get_name())

    print("Missing domains:")
    for domain in all_domains:
        print(domain.get_name())

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
        for rank, domains in rank_to_domains.items():
            out_file = os.path.join(out_dir, f"{rank}.txt")
            write_list(domains, out_file)


def main():
    args = parse_arguments()
    engine = create_engine(f"sqlite:///{args.database}")

    with Session(engine) as session:
        split_by_clade(session, args.taxonomy, args.rank, args.out)


if __name__ == "__main__":
    main()
