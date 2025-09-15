from sqlalchemy import select
from sqlalchemy.orm import Session
from sqlalchemy.sql import func

from parasect.database.build_database import AdenylationDomain, Substrate, DomainSynonym, ProteinDomainAssociation, \
    ProteinSynonym, Protein, Taxonomy
from parasect.core.chem import is_same_molecule_fingerprint, smiles_to_fingerprint
from parasect.core.taxonomy import Rank


def get_protein_names(session: Session):
    query = select(DomainSynonym)
    protein_names = set()
    domain_synonyms = session.scalars(query)
    for synonym in domain_synonyms:
        protein_names.add('.'.join(synonym.synonym.split('.')[:-1]))

    return protein_names


def sequences_are_equivalent(seq_1: str, seq_2: str, min_overlap_ratio: float = 0.8) -> bool:
    """Return True if two sequences are equivalent based on an overlap threshold, False otherwise

    :param seq_1: sequence 1
    :type seq_1: str
    :param seq_2: sequence 2
    :type seq_2: str
    :param min_overlap_ratio: proportion of the shortest sequence of the pair that needs to match
    :type min_overlap_ratio: float
    :return: True if sequences are equivalent, False otherwise
    :rtype: bool
    """

    def sequences_overlap(a: str, b: str, threshold: float) -> bool:
        """
        Return True if sequences overlap by at least a certain threshold, False otherwise

        :param a: sequence 1
        :type a: str
        :param b: sequence 2
        :type b: str
        :param threshold: proportion of shortest sequence that has to match exactly
        :type threshold: float
        :return: True if sequences overlap, False otherwise
        :rtype: bool
        """
        max_overlap = min(len(a), len(b))
        min_overlap = int(max_overlap * threshold)
        for length in range(max_overlap, min_overlap - 1, -1):
            if a[-length:] == b[:length]:
                return True
        return False

    if seq_1 in seq_2:
        return True
    if seq_2 in seq_1:
        return True

    if sequences_overlap(seq_1, seq_2, min_overlap_ratio):
        return True
    if sequences_overlap(seq_2, seq_1, min_overlap_ratio):
        return True

    return False


def get_all_substrates(session: Session) -> list[Substrate]:
    query = select(Substrate)
    substrates = []
    for substrate in session.scalars(query):
        substrates.append(substrate)

    return substrates


def get_substrates_from_smiles(session: Session, smiles: str) -> list[Substrate]:
    """

    :param session: database session
    :type session: Session
    :param smiles: SMILES string of query molecule
    :type smiles: str
    :return: list of matching substrates
    :rtype:
    """
    query = select(Substrate)
    substrates = []
    fingerprint = smiles_to_fingerprint(smiles)
    for substrate in session.scalars(query):
        if is_same_molecule_fingerprint(set(substrate.fingerprint), set(fingerprint)):
            substrates.append(substrate)

    return substrates


def get_domains_from_taxonomic_rank(session: Session, rank: Rank, rank_name: str) -> list[AdenylationDomain]:

    taxonomy_column = Rank.get_column(rank)
    query = (
        select(AdenylationDomain)
            .join(AdenylationDomain.proteins)
            .join(ProteinDomainAssociation.protein)
            .join(Protein.taxonomy)
            .where(taxonomy_column == rank_name)
            .distinct()
    )
    return list(session.scalars(query).all())


def get_substrates_from_name(session: Session, substrate_name: str) -> list[Substrate]:
    query = select(Substrate).where(func.lower(Substrate.name) == substrate_name.lower())
    return list(session.scalars(query))


def get_domains_from_synonym(session: Session, synonym: str) -> list[AdenylationDomain]:
    query = select(AdenylationDomain).join(DomainSynonym).where(DomainSynonym.synonym == synonym)
    return list(session.scalars(query))


def get_proteins_from_synonym(session: Session, synonym: str) -> list[Protein]:
    query = select(Protein).join(ProteinSynonym).where(ProteinSynonym.synonym == synonym)
    return list(session.scalars(query))


def get_proteins_from_domain_synonym(session: Session, synonym: str) -> list[Protein]:
    query = (
        select(Protein)
        .join(Protein.domains)                   # Protein → ProteinDomainAssociation
        .join(ProteinDomainAssociation.domain)   # ProteinDomainAssociation → AdenylationDomain
        .join(AdenylationDomain.synonyms)        # AdenylationDomain → DomainSynonym
        .where(DomainSynonym.synonym == synonym)
    )
    return list(session.scalars(query).unique())


def get_domains_from_protein_synonym(session: Session, synonym: str) -> list[AdenylationDomain]:
    query = (
        select(AdenylationDomain)
        .join(AdenylationDomain.proteins)        # AdenylationDomain → ProteinDomainAssociation
        .join(ProteinDomainAssociation.protein)  # ProteinDomainAssociation → Protein
        .join(Protein.synonyms)                  # Protein → ProteinSynonym
        .where(ProteinSynonym.synonym == synonym)
    )
    return list(session.scalars(query).unique())


def get_domains_from_sequence(session: Session, sequence: str) -> list[AdenylationDomain]:
    query = select(AdenylationDomain)
    domains = []
    for domain in session.scalars(query):

        if sequences_are_equivalent(sequence, domain.sequence):
            domains.append(domain)

    return domains
