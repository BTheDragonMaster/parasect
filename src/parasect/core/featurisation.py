import os
import subprocess
import warnings
from typing import List, Union

from pikachu.fingerprinting.ecfp_4 import ECFP
from pikachu.general import Structure, read_smiles

from parasect.core.abstractions import AdenylationDomain
from parasect.core.constants import PROPERTIES, HMM2_FILE
from parasect.core.parsing import (
    parse_morgan_fingerprint_file, 
    read_fasta_file,
    parse_protein_sequences_from_genbank_file,
)
from parasect.core.hmmer import (
    rename_sequences,
    reverse_renaming,
    run_hmmpfam2,
    parse_hmm2_results,
)

def bitvector_from_smiles(smiles_string: str, bitvector_file: str) -> list[int]:
    """
    Return a bitvector from a SMILES string

    Parameters
    ----------
    smiles_string: str, must be a valid SMILES string
    bitvector_file: path to file containing bitvectors on which paras and parasect have been trained

    Returns
    -------
    vector: list of int, with each int either 0 or 1, 0 denoting absence and 1 denoting prescence of bitvector. Vector
        is given in the same order as in bitvector_file

    """
    structure: Structure = read_smiles(smiles_string)

    ecfp: ECFP = ECFP(structure)
    fingerprint: set[int] = ecfp.fingerprint

    with open(bitvector_file, "r") as bitvectors:
        substructure_hashes = list(map(int, bitvectors.readline().split("\t")[2:]))

    vector: list[int] = []

    for substructure_hash in substructure_hashes:
        if substructure_hash in fingerprint:
            vector.append(1)
        else:
            vector.append(0)

    return vector


def bitvectors_from_substrate_names(substrate_names, fingerprint_file):
    """
    Return substrate names and fingerprints for all substrate names for which a fingerprint could be found

    Parameters
    ----------
    substrate_names: list of [str, ->], with each str a substrate name. If the substrate name is not in the fingerprint
        file, a warning will be printed.
    fingerprint_file: str, path to file containing precomputed fingerprints, with one substrate per row and one
        substructure per column

    Returns
    -------
    substrates: list of [str, ->], substrate names which were present in the fingerprint file
    fingerprints: list of [[int, ->], ->], with each list of integers a fingerprint. The index of each fingerprint
        matches the index of the corresponding substrate name in substrates.

    """
    substrate_to_fingerprint: dict[str, list[int]] = parse_morgan_fingerprint_file(fingerprint_file)
    substrates: list[str] = []
    fingerprints: list[list[int]] = []
    for substrate in substrate_names:
        if substrate in substrate_to_fingerprint:
            substrates.append(substrate)
            fingerprints.append(substrate_to_fingerprint[substrate])
        else:
            warnings.warn(
                f"Could not find a fingerprint for substrate name {substrate}. Excluded from analysis."
            )
    return substrates, fingerprints

def domains_to_features(
    domains: list["AdenylationDomain"], signature: str = "extended", one_hot: bool = False
):
    """
    Return a list of sequence ids and a list of corresponding feature vectors from a list of AdenylationDomain instances

    Parameters
    ----------
    domains: list of [domain, ->], with each domain an AdenylationDomain instance
    signature: str, type of signature used for featurisation. Must be 'extended' or 'short'
    one_hot: bool, use one-hot encoding if True, the 15 NRPSPredictor features otherwise

    Returns
    -------
    sequence_ids: list of [sequence_id, ->], with each sequence_id str
    feature_vectors: list of [[feature, ->], ->], with each feature float or int. Index of feature list corresponds to
        index of corresponding sequence_id in sequence_ids

    """
    sequence_ids: list[str] = []
    feature_vectors: list[Union[list[float], list[int]]] = []

    for domain in domains:
        sequence_ids.append(domain.domain_id)

        if signature == "extended":
            if not one_hot:
                feature_vector: list[float] = get_sequence_features(domain.extended_signature)
            else:
                feature_vector: list[int] = one_hot_encoding(domain.extended_signature)

        elif signature == "short":
            if not one_hot:
                feature_vector: list[float] = get_sequence_features(domain.signature)
            else:
                feature_vector: list[int] = one_hot_encoding(domain.signature)
        else:
            raise ValueError(f"Expected 'extended' or 'short'. Got {signature}.")

        feature_vectors.append(feature_vector)

    return sequence_ids, feature_vectors


def get_sequence_features(sequence: str) -> list[float]:
    """
    Return a feature list of NRPSPredictor features from a sequence

    Parameters
    ----------
    sequence: str, amino acid sequence

    Returns
    -------
    features: list of [feature, ->], with each feature float

    """
    features: list[float] = []

    for aa in sequence:
        properties = PROPERTIES[aa]
        features.extend(properties)

    return features


def one_hot_encoding(sequence: str) -> list[int]:
    """
    Return one-hot encoding of a sequence

    Parameters
    ----------
    sequence: str, amino acid sequence

    Returns
    -------
    feature_vector: list of [feature, ->], with each feature int

    """
    one_hot_categories = [
        "A",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y",
    ]

    feature_vector: list[int] = []

    for res in sequence:
        assert res in one_hot_categories or res == "-"
        for amino_acid in one_hot_categories:

            if res == amino_acid:
                feature_vector.append(1)
            else:
                feature_vector.append(0)

    return feature_vector



def hits_to_domains(id_to_hit, fasta_file, temp_dir, profile=False, verbose=False):
    hits_by_seq_id = {}
    if verbose:
        print("\tSorting hits by sequence..")

    for hit_key, hit in id_to_hit.items():
        seq_id, hit_id, hit_start, hit_end = parse_domain_id(hit_key)
        if seq_id not in hits_by_seq_id:
            hits_by_seq_id[seq_id] = []

        hits_by_seq_id[seq_id].append((hit_id, hit_start, hit_end, hit_key))

    if verbose:
        print("\tExtracting domain signatures..")

    counter = 0

    seq_id_to_domains = {}

    for seq_id, hits in hits_by_seq_id.items():
        counter += 1
        for hit_id_1, hit_start_1, hit_end_1, hit_key_1 in hits:
            if hit_id_1 == "AMP-binding":
                if seq_id not in seq_id_to_domains:
                    seq_id_to_domains[seq_id] = []
                match_found = False
                for hit_id_2, hit_start_2, hit_end_2, hit_key_2 in hits:
                    if hit_id_2 == "AMP-binding_C":
                        # Todo: check that 200 is a good cutoff score

                        if (hit_start_2 > hit_end_1) and (hit_start_2 - hit_end_1 < 200):
                            a_domain = AdenylationDomain(seq_id, hit_start_1, hit_end_2)
                            if not profile:
                                a_domain.set_domain_signatures_hmm(
                                    id_to_hit[hit_key_1], id_to_hit[hit_key_2]
                                )
                            seq_id_to_domains[seq_id].append(a_domain)
                            match_found = True
                            break

                if not match_found:
                    a_domain = AdenylationDomain(seq_id, hit_start_1, hit_end_1)
                    a_domain.set_domain_signatures_hmm(id_to_hit[hit_key_1])
                    seq_id_to_domains[seq_id].append(a_domain)

        if verbose and counter % 1000 == 0:
            print(f"\t\tProcessed {counter} proteins.")

    for seq_id, domains in seq_id_to_domains.items():
        domains.sort(key=lambda x: x.start)

    fasta = read_fasta_file(fasta_file)

    if verbose:
        print("\tSetting domain numbers..")

    for seq_id, sequence in fasta.items():
        counter = 1
        if seq_id in seq_id_to_domains:
            for a_domain in seq_id_to_domains[seq_id]:
                assert seq_id == a_domain.protein_name
                a_domain_sequence = sequence[a_domain.start : a_domain.end]
                if len(a_domain_sequence) > 100:
                    a_domain.set_sequence(a_domain_sequence)
                    a_domain.set_domain_number(counter)

                    counter += 1

    if not profile:
        if verbose:
            print("\tFiltering domains..")

        filtered_a_domains = []

        for seq_id, a_domains in seq_id_to_domains.items():
            for a_domain in a_domains:
                if (
                    a_domain.sequence
                    and a_domain.extended_signature
                    and a_domain.signature
                    and a_domain.domain_nr
                ):
                    filtered_a_domains.append(a_domain)

        if verbose:
            print("\tSorting domains by protein name..")

        filtered_a_domains.sort(key=lambda x: (x.protein_name, x.start))

    else:
        if verbose:
            print("\tExtracting domain signatures with profile alignment..")

        filtered_a_domains = []
        for seq_id, a_domains in seq_id_to_domains.items():
            for a_domain in a_domains:
                a_domain.set_domain_signatures_profile(temp_dir)
                if (
                    a_domain.sequence
                    and a_domain.extended_signature
                    and a_domain.signature
                    and a_domain.domain_nr
                ):
                    filtered_a_domains.append(a_domain)

        filtered_a_domains.sort(key=lambda x: (x.protein_name, x.start))

    return filtered_a_domains


def parse_domain_id(fasta_id):
    """
    Return id, domain type and domain location from common id

    Input:
    fasta_id: str

    Output:
    id: str, sequence id
    hit_id: str, domain id
    hit_start: int, start position of domain in protein
    hit_end: int, end position of domain in protei
    """
    seq_id, hit_id, hit_location = fasta_id.split("|")
    hit_start, hit_end = hit_location.split("-")
    hit_start = int(hit_start)
    hit_end = int(hit_end)
    return seq_id, hit_id, hit_start, hit_end


def domains_from_fasta(fasta_in, temp_dir, job_name="paras_run", profile=False, verbose=False):
    """
    Extract adomains from fasta file and write them to out_dir

    Input:
    fasta_in: str, file location of .fasta file containing aa sequences
    job_name: str, name of job

    Output:
    fasta_out: str, file location of .fasta file containing detected adomains

    """
    hmm_out = os.path.join(temp_dir, f"{job_name}.hmm_result")

    if verbose:
        print("Running HMM..")
    run_hmmpfam2(HMM2_FILE, fasta_in, hmm_out)
    if verbose:
        print("Parsing hits..")
    id_to_hit = parse_hmm2_results(hmm_out)

    if profile:
        if verbose:
            print("Processing hits (profile alignment-based active site extraction)..")
        a_domains = hits_to_domains(id_to_hit, fasta_in, temp_dir, profile=True, verbose=verbose)
    else:
        if verbose:
            print("Processing hits (hmm-based active site extraction)..")
        a_domains = hits_to_domains(id_to_hit, fasta_in, temp_dir, verbose=verbose)

    return a_domains


def get_domains(
    input_file,
    extraction_method,
    job_name,
    separator_1,
    separator_2,
    separator_3,
    verbose,
    file_type,
    temp_dir,
):
    assert extraction_method in [
        "hmm",
        "profile",
    ], f"Only supported extraction methods are 'hmm' or 'profile'. Got {extraction_method}."
    assert file_type in [
        "fasta",
        "gbk",
    ], f"Only supported file types are 'fasta' or 'gbk'. Got {file_type}."

    if file_type == "gbk":
        original_fasta = os.path.join(temp_dir, "proteins_from_genbank.fasta")
        parse_protein_sequences_from_genbank_file(input_file, original_fasta)  # TODO: ???
        mapping_file, renamed_fasta_file = rename_sequences(original_fasta, temp_dir)
    else:
        mapping_file, renamed_fasta_file = rename_sequences(input_file, temp_dir)

    if extraction_method == "profile":
        a_domains = domains_from_fasta(
            renamed_fasta_file, temp_dir=temp_dir, job_name=job_name, profile=True, verbose=verbose
        )
    elif extraction_method == "hmm":
        a_domains = domains_from_fasta(
            renamed_fasta_file, temp_dir=temp_dir, job_name=job_name, verbose=verbose
        )
    else:
        raise ValueError(
            f"Only supported extraction methods are 'hmm' or 'profile'. Got {extraction_method}."
        )

    reverse_renaming(a_domains, mapping_file)

    for a_domain in a_domains:
        a_domain.set_domain_id_separators(separator_1, separator_2, separator_3)

    return a_domains