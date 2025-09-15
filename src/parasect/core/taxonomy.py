from typing import Optional
from enum import Enum

from parasect.database.build_database import Taxonomy


class Rank(Enum):
    DOMAIN = 1
    KINGDOM = 2
    PHYLUM = 3
    CLASS = 4
    ORDER = 5
    FAMILY = 6
    GENUS = 7
    SPECIES = 8
    STRAIN = 9

    @staticmethod
    def get_column(rank: "Rank"):
        rank_to_column = {
            Rank.DOMAIN: Taxonomy.domain,
            Rank.KINGDOM: Taxonomy.kingdom,
            Rank.PHYLUM: Taxonomy.phylum,
            Rank.CLASS: Taxonomy.cls,
            Rank.ORDER: Taxonomy.order,
            Rank.FAMILY: Taxonomy.family,
            Rank.GENUS: Taxonomy.genus,
            Rank.SPECIES: Taxonomy.species,
            Rank.STRAIN: Taxonomy.strain,
        }
        return rank_to_column[rank]

    @staticmethod
    def get_rank_from_taxonomy(rank: "Rank", taxonomy: Taxonomy) -> Optional[str]:
        """
        Return specific rank from full taxonomy
        :param rank: taxonomic rank
        :type rank: Rank
        :param taxonomy: taxonomy from database
        :type taxonomy: Taxonomy
        :return: string representing taxonomic rank, None if not in database
        :rtype: str or None
        """
        rank_to_value = {Rank.DOMAIN: taxonomy.domain,
                         Rank.KINGDOM: taxonomy.kingdom,
                         Rank.PHYLUM: taxonomy.phylum,
                         Rank.CLASS: taxonomy.cls,
                         Rank.ORDER: taxonomy.order,
                         Rank.FAMILY: taxonomy.family,
                         Rank.GENUS: taxonomy.genus,
                         Rank.SPECIES: taxonomy.species,
                         Rank.STRAIN: taxonomy.strain}

        return rank_to_value[rank]

    @staticmethod
    def get_rank_type_from_string(string: str) -> Optional["Rank"]:
        """
        Convert a string to a Rank enum value.

        :param string: Taxonomic rank name (case-insensitive)
        :type string: str
        :return: Corresponding Rank, or None if invalid
        :rtype: Optional[Rank]
        """
        string_to_rank_type = {r.name: r for r in Rank}
        return string_to_rank_type.get(string.upper())
