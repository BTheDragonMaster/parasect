from sqlalchemy import select
from sqlalchemy.orm import Session
from sqlalchemy.sql import func
from parasect.database.build_database import AdenylationDomain, Substrate, DomainSynonym
from parasect.core.chem import is_same_molecule


def get_protein_names(session: Session):
    query = select(DomainSynonym)
    protein_names = set()
    domain_synonyms = session.scalars(query)
    for synonym in domain_synonyms:
        protein_names.add('.'.join(synonym.synonym.split('.')[:-1]))

    return protein_names


def sequences_are_equivalent(seq1: str, seq2: str, min_overlap_ratio: float = 0.8) -> bool:
    """
    Check if two DNA sequences represent the same sequence with partial overlaps.

    Parameters:
    - seq1, seq2: DNA sequences (strings)
    - min_overlap_ratio: minimum fraction of the shorter sequence length to consider as valid overlap

    Returns:
    - True if sequences overlap sufficiently and match in the overlapping region.
    - False otherwise.
    """

    def overlap_length(a, b):
        max_overlap = min(len(a), len(b))
        for length in range(max_overlap, 0, -1):
            if a[-length:] == b[:length]:
                return length
        return 0

    overlap_1_2 = overlap_length(seq1, seq2)
    overlap_2_1 = overlap_length(seq2, seq1)

    min_len = min(len(seq1), len(seq2))
    threshold = int(min_len * min_overlap_ratio)

    if overlap_1_2 >= threshold:
        return True
    if overlap_2_1 >= threshold:
        return True
    if seq1 in seq2:
        return True
    if seq2 in seq1:
        return True

    return False


def get_all_substrates(session: Session) -> list[Substrate]:
    query = select(Substrate)
    substrates = []
    for substrate in session.scalars(query):
        substrates.append(substrate)

    return substrates


def get_substrates_from_smiles(session: Session, smiles: str) -> list[Substrate]:
    query = select(Substrate)
    substrates = []
    for substrate in session.scalars(query):
        if is_same_molecule(substrate.smiles, smiles):
            substrates.append(substrate)

    return substrates


def get_substrates_from_name(session: Session, substrate_name: str) -> list[Substrate]:
    query = select(Substrate).where(func.lower(Substrate.name) == substrate_name.lower())
    return list(session.scalars(query))


def get_domains_from_synonym(session: Session, synonym: str, domain_index: int) -> list[AdenylationDomain]:
    query = select(AdenylationDomain).join(DomainSynonym).where(f"{DomainSynonym.synonym}.A{domain_index}" == synonym)
    return list(session.scalars(query))


def get_domains_from_sequence(session: Session, sequence: str) -> list[AdenylationDomain]:
    query = select(AdenylationDomain)
    domains = []
    for domain in session.scalars(query):

        if sequences_are_equivalent(sequence, domain.sequence):
            domains.append(domain)

    return domains
