import os
from argparse import ArgumentParser, Namespace
from collections import defaultdict
from typing import Optional

from pikachu.general import read_smiles
from sqlalchemy import create_engine
from sqlalchemy.orm import Session

from parasect.core.parsing import iterate_over_dir, parse_fasta_file, write_fasta_file, SubstrateData, \
    parse_smiles_mapping
from parasect.core.tabular import Tabular, write_tabular
from parasect.core.chem import is_same_molecule
from parasect.database.query_database import get_substrates_from_smiles, get_substrates_from_name, \
    sequences_are_equivalent, get_domains_from_sequence, get_domains_from_synonym


def parse_arguments() -> Namespace:
    """Parse arguments from command line

    :return: User-defined arguments
    :rtype: Namespace
    """
    parser = ArgumentParser(description="Add pending entries to database")
    parser.add_argument('-u', type=str, required=True, help="Path to user submissions folder")
    parser.add_argument('-db', type=str, required=True, help="Path to PARASECT database")
    parser.add_argument('-o', type=str, required=True, help="Output directory")
    parser.add_argument('-c', action="store_true", help="If given, perform checks against database")

    args = parser.parse_args()
    return args


def deduplicate_substrates(substrates: list[SubstrateData]) -> list[SubstrateData]:
    """
    :param substrates: list of substrates, each comprising substrate name and substrate SMILES
    :type substrates: list[SubstrateData]
    :return: deduplicated list of substrates
    :rtype: list[SubstrateData]
    """

    name_to_smiles = {}

    for substrate in substrates:
        if substrate.name not in name_to_smiles:
            name_to_smiles[substrate.name] = substrate.smiles
            try:
                read_smiles(substrate.smiles)
            except Exception:
                raise ValueError(f"Incorrect SMILES for substrate {substrate.name}: {substrate.smiles}. Fix before updating database")
        else:
            if not is_same_molecule(name_to_smiles[substrate.name], substrate.smiles):
                raise ValueError(f"User-submitted substrate name {substrate.name} represents two different molecules: {substrate.smiles}, {name_to_smiles[substrate.name]}. Fix before updating database.")

    names = list(name_to_smiles.keys())
    names.sort()

    for i, name_1 in enumerate(names):
        smiles_1 = name_to_smiles[name_1]
        for name_2 in names[i + 1:]:
            smiles_2 = name_to_smiles[name_2]
            if is_same_molecule(smiles_1, smiles_2):
                raise ValueError(
                    f"User-submitted molecule {smiles_1} represented by two different names: {name_1}, {name_2}. Fix before updating database.")

    deduplicated_substrates: list[SubstrateData] = []
    for name, smiles in name_to_smiles.items():
        deduplicated_substrates.append(SubstrateData(name, smiles))

    return deduplicated_substrates


def check_new_substrates(substrates: list[SubstrateData], session: Session) -> None:
    """Check new substrates against database

    :param substrates: list of new substrates
    :type substrates: list[SubstrateData]
    :param session: PARASECT database session
    :type session: Session
    """
    existing_substrates = {}
    existing_names = {}
    for substrate in substrates:
        print(f"Checking substrate {substrate.name}...")
        matching_substrates = get_substrates_from_smiles(session, substrate.smiles)
        matching_names = get_substrates_from_name(session, substrate.name)
        if matching_substrates:
            existing_substrates[substrate.name] = [m.name for m in matching_substrates]

        if matching_names:
            existing_names[substrate.name] = [m.name for m in matching_names]

    if existing_substrates:
        print("The following new substrates may already exist in the dataset under a different name:")
        for new_substrate, old_substrates in existing_substrates.items():
            print(f"{new_substrate}\t{'|'.join(old_substrates)}")

    if existing_names:
        print("The following new substrate names already exist in the dataset.")
        for new_substrate, old_substrates in existing_names.items():
            print(f"{new_substrate}\t{'|'.join(old_substrates)}")

        raise ValueError("Please remove all existing substrate names from smiles.tsv, or update names.")


def check_new_sequences(id_to_seq: dict[str, str], session: Session) -> None:
    """

    :param id_to_seq: dictionary of sequence identifiers and sequences
    :type id_to_seq: dict[str, str]
    :param session: PARASECT database session
    :type session: Session

    """
    for seq_id, seq in id_to_seq.items():
        domain_synonyms = seq_id.split('|')
        for synonym in domain_synonyms:
            if get_domains_from_synonym(session, synonym):
                raise ValueError(f"Invalid sequence id: {seq_id}: one of the domain synonyms ({synonym}) already exists in the database. Please change before updating database.")
        if get_domains_from_sequence(session, seq):
            raise ValueError(
                f"Sequence already exists in database: {seq_id}. Please remove entry before updating database.")


def deduplicate_sequences(id_to_seq: dict[str, str],
                          min_overlap_ratio: float = 0.8) -> tuple[dict[str, str], dict[str, str]]:
    """
    Remove duplicate sequences, updating the sequence id

    :param id_to_seq: dictionary of sequence identifiers and sequences
    :type id_to_seq: dict[str, str]
    :param min_overlap_ratio: fraction of the shortest sequence (start and front) which has to overlap with the largest
    sequence in order for the sequences to be considered the same.
    :type min_overlap_ratio: float
    :return: dictionary of deduplicated sequences, with updated sequence IDs to reflect all original IDs; dictionary of
    old ids to new ids
    :rtype: tuple[dict[str, str], dict[str, str]]
    """
    ids = sorted(id_to_seq.keys())

    # Union–find
    parent = {i: i for i in ids}
    rank = {i: 0 for i in ids}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra == rb:
            return
        if rank[ra] < rank[rb]:
            parent[ra] = rb
        elif rank[ra] > rank[rb]:
            parent[rb] = ra
        else:
            parent[rb] = ra
            rank[ra] += 1

    # Pairwise compare
    for i, id1 in enumerate(ids):
        for id2 in ids[i + 1:]:
            if sequences_are_equivalent(id_to_seq[id1], id_to_seq[id2], min_overlap_ratio):
                union(id1, id2)

    # Build groups
    groups = defaultdict(list)
    for i in ids:
        groups[find(i)].append(i)

    # Write output: join all IDs, pick representative
    dedup: dict[str, str] = {}
    old_to_new: dict[str, str] = {}
    for root, members in groups.items():
        members_sorted = sorted(members)
        rep = max(members, key=lambda i: len(id_to_seq[i]))

        new_id = "|".join(members_sorted)
        dedup[new_id] = id_to_seq[rep]
        for member in members_sorted:
            old_to_new[member] = new_id

    return dedup, old_to_new


def replace_illegal_character(character: str, replacement_character: str, domains_file: str, signatures_file: str,
                              extended_file: str, proteins_file: str, dataset_file: str) -> None:
    """Replace illegal characters from sequence ids in files

    :param character: illegal character to replace
    :type character: str
    :param replacement_character: character to replace the illegal character with
    :type replacement_character: str
    :param domains_file: domains.fasta file path
    :type domains_file: str
    :param signatures_file: signatures.fasta file path
    :type signatures_file: str
    :param extended_file: extended_signatures.fasta file path
    :type extended_file: str
    :param proteins_file: proteins.fasta file path
    :type proteins_file: str
    :param dataset_file: parasect_data.txt file path
    :type dataset_file: str
    """
    protein_id_to_seq = parse_fasta_file(proteins_file)
    protein_name_mappings: dict[str, str] = {}
    new_protein_to_seq: dict[str, str] = {}

    for seq_id, protein_seq in protein_id_to_seq.items():
        new_seq_id = seq_id.replace(character, replacement_character)
        protein_name_mappings[seq_id] = new_seq_id
        new_protein_to_seq[new_seq_id] = protein_seq

    write_fasta_file(new_protein_to_seq, proteins_file)

    def get_domain_id(old_id: str) -> str:
        """ Get new domain ID from old domain ID
        :param old_id: old domain ID
        :type old_id: str
        :return: new domain ID
        :rtype: str
        """
        protein_id = '.'.join(old_id.split('.')[:-1])
        domain_number = old_id.split('.A')[-1]
        new_protein_id = protein_name_mappings[protein_id]
        new_domain_id = f"{new_protein_id}.A{domain_number}"
        return new_domain_id

    def change_domain_fasta_headers(fasta_file: str) -> None:
        """
        Change domain headers based on new protein names
        :param fasta_file: path to fasta file. File will be overwritten
        :type fasta_file: str
        """
        id_to_seq = parse_fasta_file(fasta_file)
        new_id_to_seq = {}
        for domain_id, seq in id_to_seq.items():
            new_id = get_domain_id(domain_id)
            new_id_to_seq[new_id] = seq

        write_fasta_file(new_id_to_seq, fasta_file)

    change_domain_fasta_headers(domains_file)
    change_domain_fasta_headers(signatures_file)
    change_domain_fasta_headers(extended_file)

    parasect_data = Tabular(dataset_file, separator='\t')

    with open(dataset_file, 'w') as parasect_out:
        parasect_out.write("domain_id\tsequence\tspecificity\n")
        for i, domain_name in enumerate(parasect_data.rows):
            new_domain_id = get_domain_id(domain_name)
            sequence = parasect_data.get_row_value(domain_name, "sequence")
            specificity = parasect_data.get_row_value(domain_name, "specificity")
            parasect_out.write(f"{new_domain_id}\t{sequence}\t{specificity}\n")


def deduplicate_domains(domain_mappings: list[dict[str, str]],
                        parasect_datasets: list[str]) -> tuple[dict[str, str], dict[str, str], dict[str, str]]:
    """Deduplicate domains by name and by sequence

    :param domain_mappings: list of dictionaries of sequence IDs to domain sequence, with each dictionary containing
                            the domains from a single user submission
    :type domain_mappings: list[dict[str, str]
    :param parasect_datasets: list of paths to parasect data files
    :type parasect_datasets: list[str]
    :return: dictionary of sequence ID to domain sequence, dictionary of old ids to new ids, dictionary of new ids to
    specificity
    :rtype: tuple[dict[str, str], dict[str, str], dict[str, str]]
    """

    id_to_seq: dict[str, str] = {}
    id_to_spec: dict[str, str] = {}

    annotated = 0

    for i, domain_mapping in enumerate(domain_mappings):
        parasect_data = Tabular(parasect_datasets[i], separator='\t')
        for seq_id, seq in domain_mapping.items():
            specificity = parasect_data.get_row_value(seq_id, 'specificity')
            annotated += 1
            if seq_id in id_to_seq:
                if seq != id_to_seq[seq_id]:
                    raise ValueError(f"User-submitted domains of the same name found with mismatching sequences: {seq_id}. Please fix before updating database")
                if id_to_spec[seq_id] != specificity:
                    raise ValueError(
                        f"Inconsistent specificities reported for domain {seq_id}: {id_to_spec[seq_id]}, {specificity}. Fix before updating database.")
            else:
                id_to_seq[seq_id] = seq
                id_to_spec[seq_id] = specificity

    deduplicated_sequences, name_mapping = deduplicate_sequences(id_to_seq)

    check_specificities(deduplicated_sequences, id_to_spec)
    id_to_spec = remap_sequences(id_to_spec, name_mapping)

    print(f"Total annotated: {annotated}")
    print(f"Total unique: {len(deduplicated_sequences)}")

    return deduplicated_sequences, name_mapping, id_to_spec


def check_specificities(deduplicated_domains: dict[str, str], specificity_data: dict[str, str]) -> None:
    """Check that all user submissions for a single domain list the same set of specificities

    :param deduplicated_domains: dictionary of domain id to sequence
    :type deduplicated_domains: dict[str, str]
    :param specificity_data: dictionary of domain id to specificity
    :type specificity_data: dict[str, str]
    """

    for domain_name in deduplicated_domains.keys():
        if '|' in domain_name:
            domain_synonyms = domain_name.split('|')
            specificity: Optional[str] = None

            for domain_synonym in domain_synonyms:
                domain_specificity = specificity_data[domain_synonym]

                if specificity is not None:
                    if domain_specificity != specificity:
                        raise ValueError(f"Inconsistent specificities reported for domain {domain_name}: {domain_specificity}, {specificity}. Fix before updating database.")
                specificity = domain_specificity


def update_sequence_dict(id_to_seq: dict[str, str], fasta_file: str) -> None:
    """

    :param id_to_seq: dictionary of sequence id to A domain (extended) signature
    :type id_to_seq: dict[str, str]
    :param fasta_file: path to fasta file containing A domain (extended) signatures
    :type fasta_file: str
    :return:
    """
    id_to_seq_2 = parse_fasta_file(fasta_file)

    for seq_id, seq in id_to_seq_2.items():
        if seq_id in id_to_seq:
            if id_to_seq[seq_id] != seq:
                raise ValueError(f"Mismatching sequences found for {seq_id} in user submissions: {id_to_seq[seq_id]}, {seq}. Fix before updating database.")
        else:
            id_to_seq[seq_id] = seq


def remap_sequences(id_to_seq: dict[str, str], mapping: dict[str, str]) -> dict[str, str]:
    """
    Remap sequences to new identifiers

    :param id_to_seq: dictionary of old id to sequence
    :type id_to_seq: dict[str, str]
    :param mapping: dictionary of old id to new id
    :type mapping: dict[str, str]

    :return: dictionary of new id to sequence
    :rtype: dict[str, str]
    """

    new_id_to_seq: dict[str, str] = {}

    for old_id, seq in id_to_seq.items():
        new_id = mapping[old_id]
        if new_id in new_id_to_seq:
            if seq != new_id_to_seq[new_id]:
                raise ValueError(f"Mismatching signature sequences found for {new_id}: {seq}, {new_id_to_seq[new_id]}. Fix before updating database.")
        else:
            new_id_to_seq[new_id] = seq

    return new_id_to_seq


def collate_user_submissions(submission_dir: str, database_path: str, check_against_db: bool, out_dir: str) -> None:
    """
    :param submission_dir: Path to directory containing user submissions
    :type submission_dir: str
    :param database_path: Path to parasect database
    :type database_path: str
    :param out_dir: Path to output directory
    :type out_dir: str
    """
    engine = create_engine(f"sqlite:///{database_path}")
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    domains_out = os.path.join(out_dir, "domains.fasta")
    signatures_out = os.path.join(out_dir, "signatures.fasta")
    extended_out = os.path.join(out_dir, "extended_signatures.fasta")
    proteins_out = os.path.join(out_dir, "proteins.fasta")
    parasect_out = os.path.join(out_dir, "parasect_data.txt")
    smiles_out = os.path.join(out_dir, "smiles.tsv")

    all_substrates = []
    domain_mappings = []
    parasect_datasets = []
    id_to_sig: dict[str, str] = {}
    id_to_ext: dict[str, str] = {}
    protein_to_seq: dict[str, str] = {}

    for _, user_submission_path in iterate_over_dir(submission_dir, get_dirs=True):
        smiles_path = os.path.join(user_submission_path, "smiles.tsv")
        for submission_type, submission in iterate_over_dir(user_submission_path, get_dirs=True):
            domains_fasta = os.path.join(submission, "domains.fasta")
            signatures_fasta = os.path.join(submission, "signatures.fasta")
            extended_fasta = os.path.join(submission, "extended_signatures.fasta")
            parasect_data = os.path.join(submission, "parasect_data.txt")
            proteins_fasta = os.path.join(submission, "proteins.fasta")

            replace_illegal_character('|', '_', domains_fasta, signatures_fasta, extended_fasta, proteins_fasta,
                                      parasect_data)

            domain_mappings.append(parse_fasta_file(domains_fasta))
            parasect_datasets.append(parasect_data)
            update_sequence_dict(id_to_sig, signatures_fasta)
            update_sequence_dict(id_to_ext, extended_fasta)
            update_sequence_dict(protein_to_seq, proteins_fasta)

        if os.path.exists(smiles_path):
            substrates = parse_smiles_mapping(smiles_path)
            all_substrates.extend(substrates)

    print("Deduplicating domains...")
    id_to_seq, name_mapping, id_to_spec = deduplicate_domains(domain_mappings, parasect_datasets)
    print(f"{len(id_to_seq)} unique domains found.")

    print("Deduplicating substrates...")
    deduplicated_substrates = deduplicate_substrates(all_substrates)

    if check_against_db:
        with Session(engine) as session:
            print("Checking sequences against database...")
            check_new_sequences(id_to_seq, session)

            print("Checking substrates against database...")
            check_new_substrates(deduplicated_substrates, session)

    new_id_to_sig = remap_sequences(id_to_sig, name_mapping)
    new_id_to_ext = remap_sequences(id_to_ext, name_mapping)

    write_fasta_file(new_id_to_sig, signatures_out)
    write_fasta_file(new_id_to_ext, extended_out)
    write_fasta_file(id_to_seq, domains_out)
    write_fasta_file(protein_to_seq, proteins_out)
    write_tabular([id_to_seq, id_to_spec], ["domain_id", "sequence", "specificity"], parasect_out)

    with open(smiles_out, 'w') as out:
        out.write("substrate\tsmiles\n")
        for substrate in deduplicated_substrates:
            out.write(f"{substrate.name}\t{substrate.smiles}\n")


def main() -> None:
    args = parse_arguments()
    collate_user_submissions(args.u, args.db, args.c, args.o)


if __name__ == "__main__":
    main()
