from argparse import ArgumentParser, Namespace
import os

from sqlalchemy import create_engine, select
from sqlalchemy.orm import Session

from parasect.database.query_database import get_domains_from_taxonomic_rank
from parasect.core.taxonomy import Rank
from parasect.core.writers import write_list
from parasect.database.build_database import AdenylationDomain


def parse_arguments() -> Namespace:
    """Parse arguments from command line

    :return: Arguments
    :rtype: Namespace
    """
    parser = ArgumentParser(description="Retrieve taxonomy for database domains")

    parser.add_argument("-db", "--database", type=str, required=True,
                        help="Path to PARASECT database")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Path to output directory")

    args = parser.parse_args()

    return args


if __name__ == "__main__":
    args = parse_arguments()
    engine = create_engine(f"sqlite:///{args.database}")
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    with Session(engine) as session:
        bacterial_domains = [d.get_name() for d in get_domains_from_taxonomic_rank(session, Rank.DOMAIN, "Bacteria")]
        fungal_domains = [d.get_name() for d in get_domains_from_taxonomic_rank(session, Rank.KINGDOM, "Fungi")]

        all_domains = set([d.get_name() for d in list(session.scalars(select(AdenylationDomain)).all())])
        unknown = all_domains - set(bacterial_domains) - set(fungal_domains)

        bacterial_out = os.path.join(args.output, "bacterial_domains.txt")
        fungal_out = os.path.join(args.output, "fungal_domains.txt")
        unknown_out = os.path.join(args.output, "unknown.txt")

        write_list(bacterial_domains, bacterial_out)
        write_list(fungal_domains, fungal_out)
        write_list(list(unknown), unknown_out)
