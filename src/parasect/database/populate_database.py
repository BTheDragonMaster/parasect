from argparse import ArgumentParser, Namespace
import re
from typing import Any

from sqlalchemy import create_engine, select
from sqlalchemy.orm import Session

from parasect.core.parsing import load_parasect_data, parse_fasta_file, AdenylationDomainData, SubstrateData
from parasect.database.build_database import AdenylationDomain, Substrate, DomainSynonym, ProteinSynonym, Protein, \
    ProteinDomainAssociation


def parse_args() -> Namespace:
    parser = ArgumentParser(description="Populate empty SQL database with data")
    parser.add_argument("--parasect", required=True, type=str, help="Path to parasect dataset file")
    parser.add_argument("--smiles", required=True, type=str, help="Path to parasect SMILES file")
    parser.add_argument("--database", required=True, type=str, help="Path to parasect database")
    parser.add_argument("--signature", required=True, type=str, help="Path to fasta containing A-domain signatures")
    parser.add_argument("--extended", required=True, type=str,
                        help="Path to fasta containing A-domain extended signatures")
    parser.add_argument("--protein", required=True, type=str, help="Path to fasta file containing full protein sequences")
    args = parser.parse_args()

    return args


def create_protein_entries(protein_path: str, session: Session) -> tuple[list[Protein], list[ProteinSynonym]]:
    protein_entries = []
    protein_synonyms = []

    protein_to_seq = parse_fasta_file(protein_path)

    for protein, seq in protein_to_seq.items():

        synonyms = protein.split('|')

        existing_protein = session.scalar(
            select(Protein).join(Protein.synonyms).where(ProteinSynonym.synonym.in_(synonyms)))

        if existing_protein:
            entry = existing_protein
            if entry.sequence != seq:
                raise ValueError(f"Sequence of protein {protein} does not match with protein in database: {'|'.join([s.synonym for s in entry.synonyms])}")

            existing_synonyms = {s.synonym for s in entry.synonyms}
            for synonym in synonyms:
                if synonym not in existing_synonyms:
                    protein_synonyms.append(ProteinSynonym(synonym=synonym, protein=entry))
        else:
            entry = Protein(sequence=seq)

            for synonym in synonyms:
                protein_synonyms.append(ProteinSynonym(synonym=synonym, protein=entry))

        protein_entries.append(entry)

    return protein_entries, protein_synonyms


def create_domain_entries(domains: list[AdenylationDomainData], substrates: list[SubstrateData],
                          pending: bool = False) -> \
        tuple[list[AdenylationDomain], list[DomainSynonym], list[Substrate]]:

    domain_entries = []
    substrate_entries = []
    domain_synonyms = []

    for substrate in substrates:
        substrate_entries.append(Substrate(name=substrate.name, smiles=substrate.smiles))

    for domain in domains:
        domain_entry = AdenylationDomain(sequence=domain.sequence, signature=domain.signature,
                                         extended_signature=domain.extended_signature,
                                         pending=pending)
        for synonym in domain.synonyms:
            domain_synonyms.append(DomainSynonym(synonym=synonym, domain=domain_entry))

        for substrate_entry in substrate_entries:
            if substrate_entry.name in [s.name for s in domain.substrates]:
                domain_entry.substrates.append(substrate_entry)

        domain_entries.append(domain_entry)

    return domain_entries, domain_synonyms, substrate_entries


def _build_protein_lookup(proteins: list[Protein]) -> dict[str, Protein]:
    synonym_to_protein: dict[str, Protein] = {}
    for protein in proteins:
        for synonym in protein.synonyms:
            synonym_to_protein[synonym.synonym] = protein

    return synonym_to_protein


def link_domains_and_proteins(domain_entries: list[AdenylationDomain],
                              protein_entries: list[Protein],
                              session: Session) -> list[ProteinDomainAssociation]:
    links = []

    old_protein_entries: list[Protein] = list(session.scalars(select(Protein)).all())
    protein_lookup = _build_protein_lookup(old_protein_entries + protein_entries)

    for domain_entry in domain_entries:
        for domain_synonym in domain_entry.synonyms:
            synonym = domain_synonym.synonym
            synonym_parts = synonym.split('.')
            protein_name = '.'.join(synonym_parts[:-1])
            domain_number = int(synonym_parts[-1][1:])

            if protein_name in protein_lookup:
                protein = protein_lookup[protein_name]

            else:
                raise ValueError(f"Could not find protein {protein_name} in dataset.")

            start_coords = [match.start() + 1 for match in re.finditer(domain_entry.sequence, protein.sequence)]

            if len(start_coords) == 0:
                raise ValueError(f"Domain {synonym} not found in protein {protein_name}")
            else:
                for start_coord in start_coords:
                    end_coord = start_coord + len(domain_entry.sequence)
                    links.append(ProteinDomainAssociation(protein=protein,
                                                          domain=domain_entry,
                                                          domain_number=domain_number,
                                                          start=start_coord,
                                                          end=end_coord))

    return links


def populate_db(session: Session, parasect_data_path: str, smiles_path: str, signature_path: str,
                extended_path: str, protein_path: str):
    domains, substrates = load_parasect_data(parasect_data_path, smiles_path, signature_path, extended_path,
                                             session=session)

    protein_entries, protein_synonyms = create_protein_entries(protein_path, session)
    domain_entries, domain_synonyms, substrate_entries = create_domain_entries(domains, substrates)
    protein_domain_links = link_domains_and_proteins(domain_entries, protein_entries, session)

    entries: list[Any] = []

    entries.extend(domain_entries)
    entries.extend(substrate_entries)
    entries.extend(domain_synonyms)
    entries.extend(protein_entries)
    entries.extend(protein_synonyms)
    entries.extend(protein_domain_links)

    session.add_all(entries)
    session.commit()


def main():
    args = parse_args()
    engine = create_engine(f"sqlite:///{args.database}")
    with Session(engine) as session:
        try:
            populate_db(session, args.parasect, args.smiles, args.signature, args.extended, args.protein)
        except Exception as e:
            print("Something went wrong while populating the database. Rolling back.")
            print(e)
            session.rollback()
        finally:
            session.close()


if __name__ == "__main__":
    main()
