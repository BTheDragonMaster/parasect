import os.path
from argparse import ArgumentParser, Namespace

from sqlalchemy.orm import Session
from sqlalchemy import create_engine, select

from parasect.database.build_database import AdenylationDomain, Protein, ProteinDomainAssociation, ProteinSynonym
from parasect.core.featurisation import get_domains


def parse_arguments() -> Namespace:
    """
    Parse arguments for extracting A domain signatures

    :return: arguments
    :rtype: Namespace
    """
    parser = ArgumentParser(description="Extract signatures for A-domain dataset")
    parser.add_argument('-p', required=True, type=str, help="Path to input protein .fasta")
    parser.add_argument('-t', required=True, type=str, help="Path to temp dir")
    parser.add_argument('-o', required=True, type=str, help="Path to output directory")
    parser.add_argument('-db', required=True, type=str, help="Path to database")

    args = parser.parse_args()

    return args


def extract_signatures(session: Session, protein_file: str, tmp_path: str, out_path: str) -> None:
    """Extract signatures for A-domain dataset

    :param session: database session
    :type session: Session
    :param tmp_path: path to tmp parasect path for storing intermediate files
    :type tmp_path: str
    :param protein_file: Path to full protein sequences
    :type protein_file: str
    :param out_path: path to output directory
    :type out_path: str

    """

    print("Extracting domains from protein.fasta...")
    domains = get_domains(protein_file, tmp_path, "hmm", "fasta")

    # db_domain_to_domain = {}
    db_domain_to_positions = {}
    db_domain_to_relative = {}

    found_db_domains = set()
    all_db_domains = set(session.scalars(select(AdenylationDomain)).all())

    for domain in domains:
        protein_synonyms = domain.protein_name.split('|')

        unique_protein_query = select(Protein).join(Protein.synonyms).where(ProteinSynonym.synonym.in_(protein_synonyms)).distinct()
        unique_proteins = session.execute(unique_protein_query).scalars().all()
        for protein in unique_proteins:
            query = (
                select(AdenylationDomain, ProteinDomainAssociation)
                    .join(AdenylationDomain.proteins)
                    .join(ProteinDomainAssociation.protein)
                    .where(Protein.id == protein.id).distinct()
            )

            db_domains = session.execute(query).all()

            for db_domain, association in db_domains:

                name = '|'.join([s.synonym for s in association.protein.synonyms])

                positions = [p for p in domain.extended_signature_positions if p is not None]
                sig_start = min(positions) if positions else None
                sig_end = max(positions) if positions else None

                if sig_start is not None and sig_end is not None:

                    if association.start - 1 <= sig_start and association.end >= sig_end:
                        found_db_domains.add(db_domain)

                        if db_domain not in db_domain_to_positions:
                            db_domain_to_positions[db_domain] = {name: []}
                        if name not in db_domain_to_positions[db_domain]:
                            db_domain_to_positions[db_domain][name] = []

                        if db_domain.extended_signature != domain.extended_signature:
                            print(f"Mismatching signatures for domain {db_domain.get_name()}:"
                                  f"{db_domain.extended_signature}"
                                  f"{domain.extended_signature}"
                                  f"{domain.extended_signature_positions}")
                        relative_positions = [
                            pos - (association.start - 1) if pos is not None else None
                            for pos in domain.extended_signature_positions
                        ]
                        db_domain_to_positions[db_domain][name].append(domain.extended_signature_positions)
                        if db_domain in db_domain_to_relative:
                            if db_domain_to_relative[db_domain] != relative_positions:
                                print(f"Mismatching relative positions for domain {db_domain.get_name()}:"
                                      f"{db_domain_to_relative[db_domain]}"
                                      f"{relative_positions}")
                        else:
                            db_domain_to_relative[db_domain] = relative_positions

    missing_domains = all_db_domains - found_db_domains
    print("Missing domains:")
    for domain in missing_domains:
        print(domain.get_name())
    if not os.path.exists(out_path):
        os.mkdir(out_path)

    relative_positions_file = os.path.join(out_path, 'relative_positions.txt')

    with open(relative_positions_file, 'w') as relative_out:
        relative_out.write("domain_id")
        for i in range(34):
            relative_out.write(f"\tpos_{i + 1}")
        relative_out.write('\n')

        for db_domain, positions in db_domain_to_relative.items():
            domain_id = db_domain.get_name()
            relative_out.write(domain_id)
            for pos in positions:
                relative_out.write(f"\t{pos}")
            relative_out.write("\n")


def main() -> None:
    args = parse_arguments()
    engine = create_engine(f"sqlite:///{args.db}")
    try:
        with Session(engine) as session:
            extract_signatures(session, args.p, args.t, args.o)
    finally:
        session.close()


if __name__ == "__main__":
    main()
