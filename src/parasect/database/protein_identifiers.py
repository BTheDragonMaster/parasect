import re
from enum import Enum
from typing import Optional


class IdType(Enum):
    GENPEPT = "genpept"
    REFSEQ = "refseq"
    UNIPROT = "uniprot"
    UNIPARC = "uniparc"


patterns: dict[IdType, re.Pattern] = {
    IdType.GENPEPT: re.compile(r"^[A-Z]{1,3}\d{5,6}(\.\d+)?$"),
    IdType.REFSEQ: re.compile(r"^(NP|XP|YP|WP|AP)_[0-9]+(\.\d+)?$"),
    IdType.UNIPROT: re.compile(
        r"^(?:"
        r"[OPQ][0-9][A-Z0-9]{3}[0-9]"      # classic starting with O, P, Q
        r"|[A-NR-Z][0-9][A-Z0-9]{3}[0-9]"  # classic other
        r"|A0A[0-9A-Z]{7}"                 # new 10-char
        r")$"
    ),
    IdType.UNIPARC: re.compile(r"^UPI[0-9A-F]{10,12}$"),
}


def check_id_type(identifier: str) -> Optional[IdType]:
    """Return the ID type if it matches known patterns, else None."""
    for id_type, pattern in patterns.items():
        if pattern.match(identifier):
            return id_type
    return None
