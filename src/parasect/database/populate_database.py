from argparse import ArgumentParser, Namespace
import re
import traceback
from typing import Any, Union, Optional

from sqlalchemy import create_engine, select
from sqlalchemy.orm import Session

from parasect.core.parsing import parse_fasta_file, SubstrateData, \
    parse_smiles_mapping, parse_parasect_data, TaxonomyData, parse_taxonomy_file
from parasect.database.build_database import AdenylationDomain, Substrate, DomainSynonym, ProteinSynonym, Protein, \
    ProteinDomainAssociation, Taxonomy
from parasect.database.query_database import sequences_are_equivalent
from parasect.core.chem import smiles_to_fingerprint, is_same_molecule_fingerprint
from parasect.core.tabular import Tabular


def parse_args() -> Namespace:
    """Parse arguments for populating database with new entries

    :return: command line arguments
    :rtype: Namespace
    """
    parser = ArgumentParser(description="Populate empty SQL database with data")
    parser.add_argument("--parasect", required=True, type=str, help="Path to parasect dataset file")
    parser.add_argument("--smiles", type=str, default=None, help="Path to parasect SMILES file")
    parser.add_argument("--database", required=True, type=str, help="Path to parasect database")
    parser.add_argument("--signature", required=True, type=str, help="Path to fasta containing A-domain signatures")
    parser.add_argument("--extended", required=True, type=str,
                        help="Path to fasta containing A-domain extended signatures")
    parser.add_argument("--protein", required=True, type=str,
                        help="Path to fasta file containing full protein sequences")
    parser.add_argument("--taxonomy", required=True, type=str,
                        help="Path to tab-separated file mapping proteins to taxonomy")
    args = parser.parse_args()

    return args


def create_taxonomy_entries(session: Session, protein_entries: list[Protein], taxonomy_file: str) -> None:
    """Create taxonomy entries

    :param session: database session
    :type session: Session
    :param protein_entries: list of existing protein entries
    :type protein_entries: list[Protein]
    :param taxonomy_file: path to taxonomy file
    :type taxonomy_file: str
    """
    raw_taxonomy_mapping = parse_taxonomy_file(taxonomy_file)
    taxonomy_mapping: dict[str, TaxonomyData] = {}

    for raw_key, tax_data in raw_taxonomy_mapping.items():
        # Split pipe-separated protein names
        for protein_synonym in raw_key.split('|'):
            taxonomy_mapping[protein_synonym] = tax_data

    existing_taxonomies = list(session.scalars(select(Taxonomy)).all())

    # Keyed by full taxonomy
    taxonomy_lookup: dict[TaxonomyData, Taxonomy] = {TaxonomyData(t.domain, t.kingdom, t.phylum, t.cls,
                                                                  t.order, t.family, t.genus, t.species,
                                                                  t.strain): t for t in existing_taxonomies
    }

    for protein in protein_entries:
        matched_tax: Optional[TaxonomyData] = None

        # Check all synonyms for a taxonomy mapping
        for synonym_entry in protein.synonyms:
            syn = synonym_entry.synonym
            if syn in taxonomy_mapping:
                matched_tax = taxonomy_mapping[syn]
                break

        if matched_tax is None:
            # No mapping found for any synonym; skip linking
            continue

        if matched_tax in taxonomy_lookup:
            taxonomy = taxonomy_lookup[matched_tax]
        else:
            taxonomy = Taxonomy(domain=matched_tax.domain,
                                kingdom=matched_tax.kingdom,
                                phylum=matched_tax.phylum,
                                cls=matched_tax.cls,
                                order=matched_tax.order,
                                family=matched_tax.family,
                                genus=matched_tax.genus,
                                species=matched_tax.species,
                                strain=matched_tax.strain)
            taxonomy_lookup[matched_tax] = taxonomy
            session.add(taxonomy)

        protein.taxonomy = taxonomy

def create_protein_entries(session: Session, protein_path: str) -> tuple[list[Protein], list[ProteinSynonym]]:
    """Create new protein database entries from a protein.fasta file

    :param session: database session
    :type session: Session
    :param protein_path: path to protein.fasta
    :type protein_path: str
    :return: lists of protein and protein synonym database entries
    :rtype: tuple[list[Protein], list[ProteinSynonym]]
    """
    print("Creating protein entries...")
    protein_entries: list[Protein] = []
    protein_synonyms: list[ProteinSynonym] = []

    protein_to_seq: dict[str, str] = parse_fasta_file(protein_path)

    existing_proteins: list[Protein] = list(session.scalars(select(Protein)).all())
    synonym_map = _build_synonym_lookup(existing_proteins)
    sequence_map = _build_sequence_lookup(existing_proteins)

    for protein, seq in protein_to_seq.items():

        synonyms: list[str] = protein.split('|')

        existing_protein = None
        for synonym in synonyms:
            if synonym in synonym_map:
                existing_protein = synonym_map[synonym]

        if existing_protein is None:
            if seq in sequence_map:
                existing_protein = sequence_map[seq]

        if existing_protein is not None:
            if existing_protein.sequence != seq:
                raise ValueError(f"Sequence mismatch for protein {protein}")

            existing_synonyms: set[str] = {s.synonym for s in existing_protein.synonyms}
            for synonym in synonyms:
                if synonym not in existing_synonyms:
                    protein_synonyms.append(ProteinSynonym(synonym=synonym, protein=existing_protein))
                    synonym_map[synonym] = existing_protein
        else:
            new_protein: Protein = Protein(sequence=seq)
            protein_entries.append(new_protein)
            for synonym in synonyms:

                protein_synonyms.append(ProteinSynonym(synonym=synonym, protein=new_protein))
                synonym_map[synonym] = new_protein

    return protein_entries, protein_synonyms


def create_substrate_entries(session: Session, smiles_file: str) -> list[Substrate]:
    """Create new substrate database entries from a smiles.tsv file

    :param session: database session
    :type session: Session
    :param smiles_file: path to file with SMILES mapping
    :type smiles_file: str
    :return: list of new substrates
    :rtype: list[Substrate]
    """
    print("Creating substrate entries...")

    new_substrates: list[Substrate] = []
    substrates: list[SubstrateData] = parse_smiles_mapping(smiles_file)

    existing_substrates: list[Substrate] = list(session.scalars(select(Substrate)).all())
    existing_substrate_mapping = _build_substrate_lookup(existing_substrates)

    for substrate in substrates:

        fingerprint = smiles_to_fingerprint(substrate.smiles)

        if substrate.name in existing_substrate_mapping:
            existing_fingerprint = set(existing_substrate_mapping[substrate.name].fingerprint)
            if not is_same_molecule_fingerprint(fingerprint, existing_fingerprint):
                raise ValueError(f"Mismatching molecules for substrate {substrate.name} ({substrate.smiles})")
        else:
            for substrate_2 in existing_substrates:
                if is_same_molecule_fingerprint(fingerprint, set(substrate_2.fingerprint)):
                    raise ValueError(f"Substrate {substrate.name} already exists in dataset under a different name: {substrate_2.name}")
            new_substrate = Substrate(name=substrate.name, smiles=substrate.smiles, fingerprint=list(fingerprint))
            new_substrates.append(new_substrate)
            existing_substrate_mapping[substrate.name] = new_substrate

    return new_substrates


def create_domain_entries(session: Session, parasect_path: str, signature_path: str, extended_path: str,
                          substrate_entries: list[Substrate]) -> tuple[list[AdenylationDomain], list[DomainSynonym]]:
    """Create new domain database entries

    :param session: database session
    :type session: Session
    :param parasect_path: path to parasect data file
    :type parasect_path: str
    :param signature_path: path to signature .fasta file
    :type signature_path: str
    :param extended_path: path to extended signature .fasta file
    :type extended_path: str
    :param substrate_entries: list of database substrate entries (new and old)
    :type substrate_entries: list[Substrate]
    :return: list of domain entries and list of domain synonym entries
    :rtype: tuple[list[AdenylationDomain], list[DomainSynonym]]
    """
    print("Creating domain entries...")
    domain_entries: list[AdenylationDomain] = []
    domain_synonyms: list[DomainSynonym] = []

    domain_to_seq, domain_to_spec = parse_parasect_data(parasect_path)
    domain_to_sig = parse_fasta_file(signature_path)
    domain_to_ext = parse_fasta_file(extended_path)

    domains: list[str] = sorted(domain_to_seq.keys())

    existing_domains: list[AdenylationDomain] = list(session.scalars(select(AdenylationDomain)).all())
    synonym_to_domain = _build_synonym_lookup(existing_domains)
    name_to_substrate = _build_substrate_lookup(substrate_entries)
    for i, domain_id in enumerate(domains):

        sequence = domain_to_seq[domain_id]
        substrates = domain_to_spec[domain_id]

        try:
            signature = domain_to_sig[domain_id]
        except KeyError:
            raise KeyError(f"Could not find signature for domain {domain_id}")

        try:
            extended_signature = domain_to_ext[domain_id]
        except KeyError:
            raise KeyError(f"Could not find extended signature for domain {domain_id}")

        synonyms = domain_id.split('|')

        existing_domain: Optional[AdenylationDomain] = None

        for synonym in synonyms:
            if synonym in synonym_to_domain:
                existing_domain = synonym_to_domain[synonym]

        if existing_domain is None:
            for domain in existing_domains:
                if domain.sequence == sequence:
                    existing_domain = domain
                elif sequences_are_equivalent(domain.sequence, sequence):
                    print(f"Domains with equivalent but non-identical sequences found: {domain_id}, {'|'.join([s.synonym for s in domain.synonyms])}")
                    response = input("Merge domains (y/n)? ")
                    if response.lower().startswith('y'):
                        print("Merging domains.")
                        existing_domain = domain
                    else:
                        print(f"Creating new entry for {domain_id}")

        if existing_domain is not None:

            if not sequences_are_equivalent(existing_domain.sequence, sequence):
                raise ValueError(f"Mismatching sequences for domain {domain_id}")

            if existing_domain.signature != signature:
                raise ValueError(f"Mismatching signatures for domain {domain_id}")

            if existing_domain.extended_signature != extended_signature:
                raise ValueError(f"Mismatching extended signatures for domain {domain_id}")

            existing_substrate_names: set[str] = {s.name for s in existing_domain.substrates}
            if existing_substrate_names != set(substrates):
                raise ValueError(
                    f"Mismatching substrates for domain {domain_id}.\n"
                    f"Existing: {existing_substrate_names}\n"
                    f"Parsed:   {substrates}"
                )

            existing_synonyms: set[str] = {s.synonym for s in existing_domain.synonyms}

            for synonym in synonyms:
                if synonym not in existing_synonyms:
                    domain_synonyms.append(DomainSynonym(synonym=synonym, domain=existing_domain))

        else:

            domain_entry = AdenylationDomain(sequence=sequence, signature=signature,
                                             extended_signature=extended_signature)
            for synonym in synonyms:
                domain_synonyms.append(DomainSynonym(synonym=synonym, domain=domain_entry))

            for substrate in substrates:
                try:
                    domain_entry.substrates.append(name_to_substrate[substrate])
                except KeyError:
                    raise ValueError(f"No substrate found in database for {substrate} (domain: {domain_id}. \
                    Add substrate to smiles.tsv")

            existing_domains.append(domain_entry)
            domain_entries.append(domain_entry)

        if i % 100 == 0 and i != 0:
            print(f"Processed {i} domains.")

    return domain_entries, domain_synonyms


def _build_synonym_lookup(sequence_entries: Union[list[Protein], list[AdenylationDomain]]) \
        -> dict[str, Union[Protein, AdenylationDomain]]:
    """Return lookup table of sequence synonyms to sequence entries

    :param sequence_entries: list of domain entries or list of protein entries
    :type sequence_entries: Union[list[Protein], list[AdenylationDomain]]
    :return: dictionary mapping synonyms to either domain entries of protein entries
    :rtype: dict[str, Union[Protein, AdenylationDomain]]
    """
    synonym_to_sequence_entry: dict[str, Union[Protein, AdenylationDomain]] = {}
    for sequence_entry in sequence_entries:
        for synonym in sequence_entry.synonyms:
            synonym_to_sequence_entry[synonym.synonym] = sequence_entry

    return synonym_to_sequence_entry


def _build_substrate_lookup(substrate_entries: list[Substrate]) -> dict[str, Substrate]:
    """Return dictionary of substrate names to substrate entries

    :param substrate_entries: list of substrate entries
    :type substrate_entries: list[Substrate]
    :return: dictionary mapping substrte names to substrate entries
    :rtype: dict[str, Substrate]
    """

    name_to_substrate: dict[str, Substrate] = {}

    for substrate in substrate_entries:
        name_to_substrate[substrate.name] = substrate

    return name_to_substrate


def _build_sequence_lookup(protein_entries: list[Protein]) -> dict[str, Protein]:
    """Return dictionary of sequences to protein entries

    :param protein_entries: list of protein entries
    :type protein_entries: list[Protein]
    :return: dictionary mapping protein sequence to protein database entry
    :rtype: dict[str, Protein]
    """
    sequence_to_protein: dict[str, Protein] = {}

    for protein in protein_entries:
        sequence_to_protein[protein.sequence] = protein

    return sequence_to_protein


def _parse_domain_synonym(domain_synonym: str) -> tuple[str, int]:
    synonym_parts = domain_synonym.split('.')
    protein_name = '.'.join(synonym_parts[:-1])
    domain_part = synonym_parts[-1]
    if not domain_part.startswith('A'):
        raise ValueError(f"Incorrect domain name format: {domain_synonym}. Did you use pipe characters in your domain name?")
    domain_number = int(domain_part[1:])

    return protein_name, domain_number


def link_domains_and_proteins(domain_entries: list[AdenylationDomain],
                              protein_entries: list[Protein],
                              session: Session) -> list[ProteinDomainAssociation]:
    print("Mapping domains to proteins...")
    links: list[ProteinDomainAssociation] = []

    old_protein_entries: list[Protein] = list(session.scalars(select(Protein)).all())
    protein_lookup = _build_synonym_lookup(old_protein_entries + protein_entries)

    for domain_entry in domain_entries:
        seen_links: set[tuple[str, str, int]] = set()
        multiple_matches = False
        domain_numbers = set()
        start_coords_multiple: list[int] = []
        domain_synonyms_multiple: list[str] = []

        match_found = False

        for domain_synonym in domain_entry.synonyms:
            synonym = domain_synonym.synonym
            protein_name, domain_number = _parse_domain_synonym(synonym)

            domain_numbers.add(domain_number)

            if protein_name in protein_lookup:

                match_found = True
                protein = protein_lookup[protein_name]


                start_coords = [match.start() + 1 for match in re.finditer(domain_entry.sequence, protein.sequence)]

                if len(start_coords) == 0:
                    raise ValueError(f"Domain {synonym} not found in protein {protein_name}")
                elif len(start_coords) == 1:
                    start_coord = start_coords[0]
                    end_coord = start_coord + len(domain_entry.sequence)
                    link_id = (protein.get_name(), domain_entry.get_name(), start_coord)

                    if link_id not in seen_links:
                        links.append(ProteinDomainAssociation(protein=protein,
                                                              domain=domain_entry,
                                                              domain_number=domain_number,
                                                              start=start_coord,
                                                              end=end_coord))
                        seen_links.add(link_id)
                else:
                    multiple_matches = True
                    start_coords_multiple = start_coords
                    domain_synonyms_multiple.append(domain_synonym.synonym)

        if not match_found:
            raise ValueError(f"Could not find protein in dataset for domain \
            {'|'.join([s.synonym for s in domain_entry.synonyms])}.")

        if multiple_matches:
            domain_numbers = sorted(domain_numbers)
            if len(domain_numbers) != len(start_coords_multiple):
                raise ValueError(f"Number of times duplicate domain {domain_synonyms_multiple} occurs in sequence ({len(start_coords_multiple)}) does not match number of distinct domain numbers: {domain_numbers}")
            else:
                for domain_synonym in domain_synonyms_multiple:
                    protein_name, domain_number = _parse_domain_synonym(domain_synonym)
                    if protein_name in protein_lookup:
                        protein = protein_lookup[protein_name]

                        start_coord_index = domain_numbers.index(domain_number)
                        start_coord = start_coords_multiple[start_coord_index]
                        end_coord = start_coord + len(domain_entry.sequence)
                        link_id = (protein.get_name(), domain_entry.get_name(), start_coord)
                        if link_id not in seen_links:
                            links.append(ProteinDomainAssociation(protein=protein,
                                                                  domain=domain_entry,
                                                                  domain_number=domain_number,
                                                                  start=start_coord,
                                                                  end=end_coord))
                            seen_links.add(link_id)

    return links


def populate_db(session: Session, parasect_data_path: str, smiles_path: Optional[str], signature_path: str,
                extended_path: str, protein_path: str, taxonomy_path: str):

    protein_entries, protein_synonyms = create_protein_entries(session, protein_path)

    create_taxonomy_entries(session, protein_entries, taxonomy_path)

    if smiles_path is not None:
        new_substrate_entries = create_substrate_entries(session, smiles_path)
    else:
        new_substrate_entries = []

    all_substrates: list[Substrate] = new_substrate_entries + list(session.scalars(select(Substrate)).all())

    domain_entries, domain_synonyms = create_domain_entries(session, parasect_data_path, signature_path, extended_path,
                                                            all_substrates)
    protein_domain_links = link_domains_and_proteins(domain_entries, protein_entries, session)

    entries: list[Any] = []

    entries.extend(domain_entries)
    entries.extend(new_substrate_entries)
    entries.extend(domain_synonyms)
    entries.extend(protein_entries)
    entries.extend(protein_synonyms)
    entries.extend(protein_domain_links)

    session.add_all(entries)


def main():
    args = parse_args()
    engine = create_engine(f"sqlite:///{args.database}")
    with Session(engine) as session:
        try:
            populate_db(session, args.parasect, args.smiles, args.signature, args.extended, args.protein, args.taxonomy)
            session.commit()
        except Exception as e:
            print(f"[ERROR] {type(e).__name__}: {e}")
            print("Rolling back changes.")
            traceback.print_exc()
            session.rollback()
        finally:
            session.close()


if __name__ == "__main__":
    main()
