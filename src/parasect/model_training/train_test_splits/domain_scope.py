from enum import Enum
import os
from argparse import ArgumentParser, Namespace

from sqlalchemy import select, create_engine
from sqlalchemy.orm import Session

from parasect.database.build_database import AdenylationDomain
from parasect.database.query_database import get_domains_from_taxonomic_rank
from parasect.core.taxonomy import Rank
from parasect.core.parsing import parse_list
from parasect.core.writers import write_list


def parse_arguments() -> Namespace:
    parser = ArgumentParser(description="Obtain fungal and bacterial domain lists")
    parser.add_argument('-i', '--input_domains', default=None, type=str,
                        help="If given, only do the analysis for domains in this file")
    parser.add_argument('-db', '--database', type=str, required=True, help="Path to PARASECT database")
    parser.add_argument('-o', '--output', required=True, type=str, help="Path to output directory")

    args = parser.parse_args()
    return args


class DomainScope(Enum):
    BACTERIAL_ONLY = 1
    FUNGAL_ONLY = 2
    ALL = 3

    @staticmethod
    def get_domains(session: Session, included_domains: "DomainScope"):
        if included_domains == DomainScope.ALL:
            domains = list(session.scalars(select(AdenylationDomain)).all())
        elif included_domains == DomainScope.FUNGAL_ONLY:
            domains = get_domains_from_taxonomic_rank(session, Rank.KINGDOM, "Fungi")
        elif included_domains == DomainScope.BACTERIAL_ONLY:
            domains = get_domains_from_taxonomic_rank(session, Rank.DOMAIN, "Bacteria")
        else:
            raise ValueError(f"Unknown domain scope: {DomainScope.name}")

        return domains


if __name__ == "__main__":
    args = parse_arguments()
    engine = create_engine(f"sqlite:///{args.database}")
    if not os.path.exists(args.output):
        os.mkdir(args.output)
    with Session(engine) as session:
        fungal_domain_names = [domain.get_name() for domain in
                               DomainScope.get_domains(session, DomainScope.FUNGAL_ONLY)]
        bacterial_domain_names = [domain.get_name() for domain in DomainScope.get_domains(session, DomainScope.BACTERIAL_ONLY)]
        if args.input_domains is not None:
            domain_names = parse_list(args.input_domains)

            fungal_domain_names = list(set(domain_names).intersection(set(fungal_domain_names)))
            bacterial_domain_names = list(set(domain_names).intersection(set(bacterial_domain_names)))

        write_list(bacterial_domain_names, os.path.join(args.output, "bacterial.txt"))
        write_list(fungal_domain_names, os.path.join(args.output, "fungal.txt"))
