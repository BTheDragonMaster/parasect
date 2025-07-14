# -*- coding: utf-8 -*-

"""Module for featurising adenylation domains."""

import logging
import os
from typing import Dict, List, Tuple

from Bio.SearchIO._model.hsp import HSP

from parasect.core.constants import HMM2_FILE, HMM3_FILE, PROPERTIES
from parasect.core.domain import AdenylationDomain
from parasect.core.hmmer import parse_hmm_results, rename_sequences, reverse_renaming, run_hmmpfam2, run_hmmscan
from parasect.core.parsing import parse_fasta_file, parse_genbank_file


def get_domain_features(amino_acid_sequence: str) -> List[float]:
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


def _hits_to_domains(
    id_to_hit: Dict[str, HSP],
    path_in_fasta_file: str,
    path_temp_dir: str,
    use_profile_alignment: bool = False,
        hmm_version: int = 2
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

                            if hmm_version == 2 and not use_profile_alignment:
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

                    if hmm_version == 2 and not use_profile_alignment:

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

                a_domain_sequence = sequence[a_domain.start:a_domain.end]

                if len(a_domain_sequence) > 100:
                    a_domain.set_sequence(a_domain_sequence)
                    a_domain.set_domain_number(counter)

                    counter += 1

    if hmm_version == 2 and not use_profile_alignment:
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
                if use_profile_alignment:
                    a_domain.set_domain_signatures_profile(path_temp_dir)
                    if (
                        a_domain.sequence
                        and a_domain.extended_signature
                        and a_domain.signature
                        and a_domain.domain_nr
                    ):
                        filtered_a_domains.append(a_domain)
                else:
                    if a_domain.sequence and a_domain.domain_nr:
                        filtered_a_domains.append(a_domain)

        filtered_a_domains.sort(key=lambda x: (x.protein_name, x.start))

    return filtered_a_domains


# TODO: keep domains detected by hmmer3 that were not detected by hmmer2
# TODO: enable profile alignment if hmmer3 detected a domain that hmmer2 didn't

def update_hmmer2_domain_sequences(hmmer2_domains: list[AdenylationDomain],
                                   hmmer3_domains: list[AdenylationDomain],
                                   path_in_fasta_file: str) -> None:
    """
    Update hmmer2 domains with hmmer3 domain sequences, such that the longest detected sequence is maintained.

    :param hmmer2_domains: list of [AdenylationDomain, ->], list of adenylation domains detected by HMMer2. These domains also
    contain information on signatures and extended signatures. Typically contain short sequences
    :type hmmer2_domains: list[AdenylationDomain]
    :param hmmer3_domains: list of [AdenylationDomain, ->], list of adenylation domains detected by HMMer2. These domains do
    not contain information on signatures and extended signatures. Typically contain full-length sequences. Use these
    sequences to update
    :type hmmer3_domains: list[AdenylationDomain]
    :param path_in_fasta_file: Path to fasta file.
    :type path_in_fasta_file: str
    """

    fasta = parse_fasta_file(path_in_fasta_file)

    for domain_1 in hmmer2_domains:
        for domain_2 in hmmer3_domains:
            if domain_1.protein_name == domain_2.protein_name and domain_1.domains_overlap(domain_2, threshold=50):
                sequence_altered = False
                if domain_2.end > domain_1.end:
                    domain_1.end = domain_2.end
                    sequence_altered = True
                if domain_2.start < domain_1.start:
                    domain_1.start = domain_2.start
                    sequence_altered = True

                if sequence_altered:
                    if domain_1.protein_name not in fasta:
                        raise ValueError("Mismatching protein names")
                    domain_1.set_sequence(fasta[domain_1.protein_name][domain_1.start:domain_1.end])


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
    hmm2_out = os.path.join(path_temp_dir, "run.hmm2_result")
    hmm3_out = os.path.join(path_temp_dir, "run.hmm3_result")

    run_hmmpfam2(HMM2_FILE, path_in_fasta_file, hmm2_out)
    run_hmmscan(HMM3_FILE, path_in_fasta_file, hmm3_out)

    id_to_hit_2 = parse_hmm_results(hmm2_out, 2)
    id_to_hit_3 = parse_hmm_results(hmm3_out, 3)

    if use_profile_alignment:
        # processing hits (profile alignment-based active site extraction)
        a_domains = _hits_to_domains(
            id_to_hit=id_to_hit_3,
            path_in_fasta_file=path_in_fasta_file,
            path_temp_dir=path_temp_dir,
            use_profile_alignment=True,
            hmm_version=3
        )
    else:
        # processing hits (hmm-based active site extraction)
        a_domains = _hits_to_domains(
            id_to_hit=id_to_hit_2,
            path_in_fasta_file=path_in_fasta_file,
            path_temp_dir=path_temp_dir,
            use_profile_alignment=False,
            hmm_version=2
        )
        a_domains_3 = _hits_to_domains(
            id_to_hit=id_to_hit_3,
            path_in_fasta_file=path_in_fasta_file,
            path_temp_dir=path_temp_dir,
            use_profile_alignment=False,
            hmm_version=3)

        update_hmmer2_domain_sequences(a_domains, a_domains_3, path_in_fasta_file)

    return a_domains


def get_domains(
    path_in: str,
    path_temp_dir: str,
    extraction_method: str,
    file_type: str,
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

    return a_domains
