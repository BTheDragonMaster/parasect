from enum import Enum

from sqlalchemy import select

from parasect.database.build_database import AdenylationDomain
from parasect.database.query_database import get_domains_from_taxonomic_rank
from parasect.core.taxonomy import Rank


class DomainScope(Enum):
    BACTERIAL_ONLY = 1
    FUNGAL_ONLY = 2
    ALL = 3

    @staticmethod
    def get_domains(session, included_domains: "DomainScope"):
        if included_domains == DomainScope.ALL:
            domains = list(session.scalars(select(AdenylationDomain)).all())
        elif included_domains == DomainScope.FUNGAL_ONLY:
            domains = get_domains_from_taxonomic_rank(session, Rank.KINGDOM, "Fungi")
        elif included_domains == DomainScope.BACTERIAL_ONLY:
            domains = get_domains_from_taxonomic_rank(session, Rank.DOMAIN, "Bacteria")
        else:
            raise ValueError(f"Unknown domain scope: {DomainScope.name}")

        return domains
    