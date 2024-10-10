# -*- coding: utf-8 -*-

"""Module for featurising adenylation domains."""

import logging
import os
from typing import Dict, List, Tuple

from Bio.SearchIO._model.hsp import HSP

from parasect.core.abstractions import AdenylationDomain
from parasect.core.constants import HMM2_FILE, PROPERTIES
from parasect.core.hmmer import parse_hmm2_results, rename_sequences, reverse_renaming, run_hmmpfam2
from parasect.core.parsing import parse_fasta_file, parse_genbank_file


def _get_sequence_features(amino_acid_sequence: str) -> List[float]:
    """Return a feature list of NRPSPredictor features from a sequence.

    :param amino_acid_sequence: Amino acid sequence.
    :type amino_acid_sequence: str
    :return: features. List of NRPSPredictor features.
    :rtype: List[float]
    """
    features = []

    for amino_acid_id in amino_acid_sequence:
        properties = PROPERTIES[amino_acid_id]
        features.extend(properties)

    return features


def domains_to_features(domains: List[AdenylationDomain]) -> Tuple[List[str], List[List[float]]]:
    """Featurise a list of AdenylationDomain instances.

    :param domains: AdenylationDomain instances.
    :type domains: List[AdenylationDomain]
    :return: sequence_ids, feature_vectors. Index of sequence_ids corresponds to
        index of feature_vectors.
    :rtype: Tuple[List[str], List[List[float]]]
    :raises ValueError: if domain id is not set.
    :raises ValueError: if extended domain signature is not set.
    """
    # initialise output lists
    sequence_ids = []
    feature_vectors = []

    for domain in domains:

        # append domain id to sequence_ids
        if domain.domain_id is None:
            raise ValueError("domain id not set")
        sequence_ids.append(domain.domain_id)

        if domain.extended_signature is None:
            raise ValueError("extended signature not set")
        feature_vector = _get_sequence_features(domain.extended_signature)

        # append feature vector to feature_vectors
        feature_vectors.append(feature_vector)

    return sequence_ids, feature_vectors


def _hits_to_domains(
    id_to_hit: Dict[str, HSP],
    path_in_fasta_file: str,
    path_temp_dir: str,
    use_profile_alignment: bool = False,
) -> List[AdenylationDomain]:
    """Extract adenylation domains from HMM hits.

    :param id_to_hit: Dictionary of HMM hits.
    :type id_to_hit: Dict[str, HSP]
    :param path_in_fasta_file: Path to fasta file.
    :type path_in_fasta_file: str
    :param path_temp_dir: Path to temporary directory.
    :type path_temp_dir: str
    :param use_profile_alignment: Use profile alignment if True, HMM otherwise.
    :type use_profile_alignment: bool
    :return: filtered_a_domains. List of AdenylationDomain instances.
    :rtype: List[AdenylationDomain]
    :raises ValueError: If protein name mismatch.
    """
    logger = logging.getLogger(__name__)

    logger.debug("sorting hits by sequence ...")

    hits_by_seq_id: Dict[str, List[Tuple[str, int, int, str]]] = {}
    for hit_key in id_to_hit.keys():

        # parse domain ID
        seq_id, hit_id, hit_location = hit_key.split("|")
        hit_start_str, hit_end_str = hit_location.split("-")
        hit_start = int(hit_start_str)
        hit_end = int(hit_end_str)

        if seq_id not in hits_by_seq_id:
            hits_by_seq_id[seq_id] = []

        hits_by_seq_id[seq_id].append((hit_id, hit_start, hit_end, hit_key))

    logger.debug("extracting domain signatures ...")

    counter = 0
    seq_id_to_domains: Dict[str, List[AdenylationDomain]] = {}
    for seq_id, hits in hits_by_seq_id.items():
        counter += 1

        for hit_id_1, hit_start_1, hit_end_1, hit_key_1 in hits:

            if hit_id_1 == "AMP-binding":
                if seq_id not in seq_id_to_domains:
                    seq_id_to_domains[seq_id] = []

                match_found = False
                for hit_id_2, hit_start_2, hit_end_2, hit_key_2 in hits:

                    if hit_id_2 == "AMP-binding_C":
                        if hit_start_2 > hit_end_1 and (hit_start_2 - hit_end_1) < 200:
                            a_domain = AdenylationDomain(
                                protein_name=seq_id, domain_start=hit_start_1, domain_end=hit_end_2
                            )

                            if not use_profile_alignment:
                                a_domain.set_domain_signatures_hmm(
                                    hit_n_terminal=id_to_hit[hit_key_1],
                                    hit_c_terminal=id_to_hit[hit_key_2],
                                )

                            seq_id_to_domains[seq_id].append(a_domain)
                            match_found = True
                            break

                if not match_found:
                    a_domain = AdenylationDomain(
                        protein_name=seq_id, domain_start=hit_start_1, domain_end=hit_end_1
                    )

                    a_domain.set_domain_signatures_hmm(
                        hit_n_terminal=id_to_hit[hit_key_1],
                        hit_c_terminal=None,
                    )
                    seq_id_to_domains[seq_id].append(a_domain)

        if counter % 1000 == 0:
            logger.debug(f"processed {counter} proteins ...")

    # sort domains by start position
    for domains in seq_id_to_domains.values():
        domains.sort(key=lambda x: x.start)

    fasta = parse_fasta_file(path_in_fasta_file)

    logger.debug("setting domain numbers ...")

    for seq_id, sequence in fasta.items():
        counter = 1

        if seq_id in seq_id_to_domains:
            for a_domain in seq_id_to_domains[seq_id]:

                if seq_id != a_domain.protein_name:
                    raise ValueError("Protein name mismatch")

                a_domain_sequence = sequence[a_domain.start : a_domain.end]

                if len(a_domain_sequence) > 100:
                    a_domain.set_sequence(a_domain_sequence)
                    a_domain.set_domain_number(counter)

                    counter += 1

    if not use_profile_alignment:
        logger.debug("filtering domains ...")

        filtered_a_domains = []
        for a_domains in seq_id_to_domains.values():
            for a_domain in a_domains:
                if (
                    a_domain.sequence
                    and a_domain.extended_signature
                    and a_domain.signature
                    and a_domain.domain_nr
                ):
                    filtered_a_domains.append(a_domain)

        logger.debug("sorting domains by protein name ...")

        filtered_a_domains.sort(key=lambda x: (x.protein_name, x.start))

    else:
        logger.debug("extracting domain signatures with profile alignment ...")

        filtered_a_domains = []
        for a_domains in seq_id_to_domains.values():
            for a_domain in a_domains:
                a_domain.set_domain_signatures_profile(path_temp_dir)
                if (
                    a_domain.sequence
                    and a_domain.extended_signature
                    and a_domain.signature
                    and a_domain.domain_nr
                ):
                    filtered_a_domains.append(a_domain)

        filtered_a_domains.sort(key=lambda x: (x.protein_name, x.start))

    return filtered_a_domains


def _domains_from_fasta(
    path_in_fasta_file: str,
    path_temp_dir: str,
    use_profile_alignment: bool = False,
) -> List[AdenylationDomain]:
    """Extract adomains from amino acid sequences in a fasta file.

    :param path_in_fasta_file: Path to fasta file.
    :type path_in_fasta_file: str
    :param path_temp_dir: Path to temporary directory.
    :type path_temp_dir: str
    :param use_profile_alignment: Use profile alignment if True, HMM otherwise.
    :type use_profile_alignment: bool
    :return: a_domains. List of AdenylationDomain instances.
    :rtype: List[AdenylationDomain]
    """
    # run HMM search and parse results
    hmm_out = os.path.join(path_temp_dir, "run.hmm_result")
    run_hmmpfam2(HMM2_FILE, path_in_fasta_file, hmm_out)
    id_to_hit = parse_hmm2_results(hmm_out)

    if use_profile_alignment:
        # processing hits (profile alignment-based active site extraction)
        a_domains = _hits_to_domains(
            id_to_hit=id_to_hit,
            path_in_fasta_file=path_in_fasta_file,
            path_temp_dir=path_temp_dir,
            use_profile_alignment=True,
        )
    else:
        # processing hits (hmm-based active site extraction)
        a_domains = _hits_to_domains(
            id_to_hit=id_to_hit,
            path_in_fasta_file=path_in_fasta_file,
            path_temp_dir=path_temp_dir,
            use_profile_alignment=False,
        )

    return a_domains


def get_domains(
    path_in: str,
    path_temp_dir: str,
    extraction_method: str,
    file_type: str,
    separator_1: str,
    separator_2: str,
    separator_3: str,
) -> List[AdenylationDomain]:
    """Extract adenylation domains from a fasta or genbank file.

    :param path_in: Path to input file.
    :type path_in: str
    :param path_temp_dir: Path to temporary directory.
    :type path_temp_dir: str
    :param extraction_method: Extraction method. Must be 'hmm' or 'profile'.
    :type extraction_method: str
    :param file_type: File type. Must be 'fasta' or 'gbk'.
    :type file_type: str
    :param separator_1: Separator 1.
    :type separator_1: str
    :param separator_2: Separator 2.
    :type separator_2: str
    :param separator_3: Separator 3.
    :type separator_3: str
    :return: a_domains. List of AdenylationDomain instances.
    :rtype: List[AdenylationDomain]
    :raises ValueError: If extraction_method or file_type is invalid.
    """
    # check input
    if extraction_method not in ["hmm", "profile"]:
        ValueError(
            f"Only supported extraction methods are 'hmm' or 'profile'. Got {extraction_method}."
        )  # noqa: E501

    if file_type not in ["fasta", "gbk"]:
        ValueError(
            f"Only supported file types are 'fasta' or 'gbk'. Got {file_type}."
        )  # noqa: E501

    if file_type == "gbk":
        # parse genbank file
        original_fasta = os.path.join(path_temp_dir, "proteins_from_genbank.fasta")
        parse_genbank_file(
            path_in=path_in, path_out=original_fasta
        )  # creates a fasta file at path_out
        mapping_file, renamed_fasta_file = rename_sequences(
            path_in=original_fasta, path_out=path_temp_dir
        )

    else:
        # parse fasta file
        mapping_file, renamed_fasta_file = rename_sequences(path_in=path_in, path_out=path_temp_dir)

    if extraction_method == "profile":
        a_domains = _domains_from_fasta(
            path_in_fasta_file=renamed_fasta_file,
            path_temp_dir=path_temp_dir,
            use_profile_alignment=True,
        )
    elif extraction_method == "hmm":
        a_domains = _domains_from_fasta(
            path_in_fasta_file=renamed_fasta_file,
            path_temp_dir=path_temp_dir,
            use_profile_alignment=False,
        )
    else:
        msg = f"Only supported extraction methods are 'hmm' or 'profile'. Got {extraction_method}."  # noqa: E501
        raise ValueError(msg)

    reverse_renaming(adenylation_domains=a_domains, path_in_mapping_file=mapping_file)

    for a_domain in a_domains:
        a_domain.set_domain_id_separators(separator_1, separator_2, separator_3)

    return a_domains
