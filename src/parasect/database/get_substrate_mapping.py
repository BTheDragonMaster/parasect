from argparse import ArgumentParser, Namespace
from sqlalchemy import create_engine, select
from sqlalchemy.orm import Session

from parasect.database.build_database import AdenylationDomain


def parse_arguments() -> Namespace:
    """ Parse command line arguments
    :return: arguments
    :rtype: Namespace
    """
    parser = ArgumentParser(description="Retrieve flatfiles from database")
    parser.add_argument('-db', required=True, type=str, help="Path to database")
    parser.add_argument('-o', required=True, type=str, help="Path to output file")

    args = parser.parse_args()
    return args


def write_domain_to_substrate_mapping(session: Session, out_file: str) -> None:

    with open(out_file, 'w') as out:
        out.write("domain_id\tspecificity\n")
        domains = list(session.scalars(select(AdenylationDomain)).all())
        domains.sort(key=lambda x: x.get_name())
        for domain in domains:
            substrate_string = '|'.join([s.name for s in domain.substrates])
            out.write(f"{domain.get_name()}\t{substrate_string}\n")


def main():
    args = parse_arguments()
    engine = create_engine(f"sqlite:///{args.db}")
    with Session(engine) as session:
        try:
            write_domain_to_substrate_mapping(session, args.o)
        finally:
            session.close()

if __name__ == "__main__":
    main()