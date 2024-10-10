# -*- coding: utf-8 -*-

"""Module for defining abstract classes and interfaces."""

import os
from typing import List, Optional, Tuple

from Bio.SearchIO._model.hsp import HSP

from parasect.core.constants import (
    ALIGNMENT_FILE,
    HMM2_POSITIONS_EXTENDED_SIGNATURE,
    HMM2_POSITIONS_SIGNATURE,
    POSITIONS_EXTENDED_SIGNATURE,
    POSITIONS_SIGNATURE,
)
from parasect.core.muscle import run_muscle
from parasect.core.parsing import parse_fasta_file


def _get_reference_positions(positions: List[int], aligned_reference: str) -> List[int]:
    """Adjust a list of positions to account for gaps in the reference sequence.

    :param positions: List of positions of interest in the reference sequence.
    :type positions: List[int]
    :param aligned_reference: The aligned reference sequence.
    :type aligned_reference: str
    :return: List of positions, each >= the original position.
    :rtype: List[int]
    """
    # adjust position of interest to account for gaps in the ref sequence alignment
    positions = []
    position = 0

    # iterate over the aligned reference sequence
    for amino_acid_idx, amino_acid_id in enumerate(aligned_reference):
        if amino_acid_id != "-":
            if position in positions:
                positions.append(amino_acid_idx)

            position += 1

    return positions


def _get_reference_positions_hmm(
    query_sequence: str, reference_sequence: str, reference_positions: List[int]
) -> Optional[str]:
    """Extract the given positions from a query alignment.

    :param query_sequence: The aligned query sequence.
    :type query_sequence: str
    :param reference_sequence: The aligned reference sequence.
    :type reference_sequence: str
    :param reference_positions: The positions of interest in the unaligned reference.
    :type reference_positions: List[int]
    :return: A string containing the sequence elements at the adjusted reference,
        or None if the positions could not be extracted.
    :rtype: str
    :raises ValueError: If the reference sequence is too short.

    .. note:: Positions are adjusted to account for gaps in the reference sequence.
    .. note:: This function assumes that the reference sequence is shorter
        than the query sequence.
    .. note:: From antiSMASH 7.0.0.
    """
    # check if the reference is too short
    if len(reference_sequence) < len(reference_positions):
        raise ValueError(
            f"reference sequence is too short: {len(reference_sequence)} < {len(reference_positions)}"
        )  # noqa: E501

    # check if the reference sequence is the same length as the query sequence
    # or if the reference sequence is shorter than the query sequence
    if not (
        len(reference_sequence) == len(query_sequence)
        or len(reference_sequence) < len(query_sequence)
    ):
        raise ValueError(
            f"reference sequence is too long: {len(reference_sequence)} != {len(query_sequence)}"
        )

    # adjust position of interest to account for gaps in the ref sequence alignment
    positions = []
    position_skipping_gaps = 0  # position in the reference sequence

    for amino_acid_idx, amino_acid_id in enumerate(reference_sequence):

        if amino_acid_id in "-.":
            continue

        if position_skipping_gaps in reference_positions:
            positions.append(amino_acid_idx)

        # increment the position in the reference sequence
        position_skipping_gaps += 1

    # check if the number of positions extracted is the same as the number of reference
    # positions to extract
    if len(positions) != len(reference_positions):
        return None

    # extract positions from query sequence
    return "".join([query_sequence[i] for i in positions])


def _align_adenylation_domain(
    domain_name: str,
    domain_sequence: str,
    alignment_file: str,
    path_temp_dir: str,
) -> Tuple[str, str]:
    """Align adanylation domain to database of adenylation domains.

    :param domain_name: The name of the domain.
    :type domain_name: str
    :param domain_sequence: The sequence of the domain.
    :type domain_sequence: str
    :param alignment_file: The path to the alignment file.
    :type alignment_file: str
    :param path_temp_dir: The path to the temporary directory.
    :type path_temp_dir: str
    :return: A tuple containing the aligned domain and the aligned reference domain.
    :rtype: Tuple[str, str]
    """
    # create paths for temporary files in temp directory
    temp_in = os.path.join(path_temp_dir, "temp_in_alignment.fasta")
    temp_out = os.path.join(path_temp_dir, "temp_out_alignment.fasta")

    with open(temp_in, "w") as temp:
        temp.write(f">{domain_name}\n{domain_sequence}")

    # align sequence to all adenylation domains in database with muscle
    run_muscle(path_in=temp_in, path_in_alignment=alignment_file, path_out=temp_out)
    aligned_domains = parse_fasta_file(temp_out)

    # aligned sequence of domain
    aligned_domain = aligned_domains[domain_name]

    # aligned sequence of 1AMU reference sequence
    aligned_reference = aligned_domains["BAA00406.1.A1"]

    return aligned_domain, aligned_reference


class AdenylationDomain:
    """Class for representing adenylation domains."""

    def __init__(self, protein_name: str, domain_start: int, domain_end: int) -> None:
        """Initialize an AdenylationDomain object.

        :param protein_name: The name of the protein.
        :type protein_name: str
        :param domain_start: The start position of the domain.
        :type domain_start: int
        :param domain_end: The end position of the domain.
        :type domain_end: int
        """
        self.protein_name = protein_name
        self.start = domain_start
        self.end = domain_end

        self.sequence: Optional[str] = None
        self.signature: Optional[str] = None
        self.extended_signature: Optional[str] = None
        self.domain_nr = 0
        self._domain_id: Optional[str] = None

    @property
    def domain_id(self) -> str:
        """Return the domain ID.

        :return: The domain ID.
        :rtype: str
        :raises Exception: If the domain ID is not defined.
        """
        if self._domain_id is None:
            raise Exception("domain ID not defined")

        return self._domain_id

    def set_domain_number(self, domain_nr: int) -> None:
        """Set the domain number.

        :param domain_nr: The domain number.
        :type domain_nr: int

        .. note:: This function modifies the domain number attribute.
        """
        self.domain_nr = domain_nr

    def set_sequence(self, sequence: str) -> None:
        """Set the sequence of the domain.

        :param sequence: The sequence of the domain.
        :type sequence: str

        .. note:: This function modifies the sequence attribute.
        """
        self.sequence = sequence

    def set_domain_id_separators(
        self, separator_1: str, separator_2: str, separator_3: str
    ) -> None:
        """Set the domain ID separators.

        :param separator_1: The first separator. This is used to separate the protein
            name from the domain number, and the domain number from the start and end
            positions.
        :type separator_1: str
        :param separator_2: The second separator. This is used to separate the domain
            name 'domain' from the domain number.
        :type separator_2: str
        :param separator_3: The third separator. This is used to separate the start and
            end positions.
        :type separator_3: str
        :raises Exception: If the protein name or domain number is not defined.

        .. note:: This functions modifies the domain ID attribute.
        """
        if not self.protein_name or not self.domain_nr:
            raise Exception("protein name and domain number need to be defined first")

        self._domain_id = (
            f"{self.protein_name}{separator_1}"
            f"domain{separator_2}{self.domain_nr}{separator_1}"
            f"{self.start}{separator_3}{self.end}"
        )

    def set_domain_signatures_hmm(
        self, hit_n_terminal: HSP, hit_c_terminal: Optional[HSP] = None
    ) -> None:
        """Extract (extended) signatures from adenylation domains using HMM profile.

        :param hit_n_terminal: The hit object for the N-terminal domain.
        :type hit_n_terminal: HSP
        :param hit_c_terminal: The hit object for the C-terminal domain.
        :type hit_c_terminal: Optional[HSP]

        .. note:: This function modifies the signature and extended signature attributes.
        """
        valid = {
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
            "-",
        }

        signature_positions = HMM2_POSITIONS_SIGNATURE
        extended_signature_positions = HMM2_POSITIONS_EXTENDED_SIGNATURE
        position_k = [36]  # hmm2 position k

        profile = hit_n_terminal.aln[1].seq
        query = hit_n_terminal.aln[0].seq
        offset = hit_n_terminal.hit_start

        signature = _get_reference_positions_hmm(
            query_sequence=query,
            reference_sequence=profile,
            reference_positions=[p - offset for p in signature_positions],
        )

        if signature and all([char in valid for char in signature]):
            self.signature = signature

        lysine = None
        if hit_c_terminal:
            profile_c = hit_c_terminal.aln[1].seq
            query_c = hit_c_terminal.aln[0].seq
            offset_c = hit_c_terminal.hit_start

            lysine = _get_reference_positions_hmm(
                query_sequence=query_c,
                reference_sequence=profile_c,
                reference_positions=[p - offset_c for p in position_k],
            )

        if self.signature:
            if lysine and lysine in valid and lysine != "-":
                self.signature += lysine
            else:
                self.signature += "K"

        extended_signature = _get_reference_positions_hmm(
            query_sequence=query,
            reference_sequence=profile,
            reference_positions=[p - offset for p in extended_signature_positions],
        )

        if extended_signature and all([char in valid for char in extended_signature]):
            self.extended_signature = extended_signature

    def set_domain_signatures_profile(self, path_temp_dir: str) -> None:
        """Extract (extended) signatures from adenylation domains using profile alignment.

        :param path_temp_dir: The path to the temporary directory.
        :type path_temp_dir: str
        :raises Exception: If the sequence is not defined.

        .. note:: This function modifies the signature and extended signature attributes.
        """
        if self.sequence is None:
            raise Exception("sequence needs to be defined first")

        aligned_domain, aligned_reference = _align_adenylation_domain(
            domain_name="DOMAIN_TO_QUERY",
            domain_sequence=self.sequence,
            alignment_file=ALIGNMENT_FILE,
            path_temp_dir=path_temp_dir,
        )

        aligned_positions_signature = _get_reference_positions(
            positions=POSITIONS_SIGNATURE, aligned_reference=aligned_reference
        )

        aligned_positions_extended_signature = _get_reference_positions(
            positions=POSITIONS_EXTENDED_SIGNATURE, aligned_reference=aligned_reference
        )

        signature = []
        for position in aligned_positions_signature:
            signature.append(aligned_domain[position])
        self.signature = "".join(signature)

        extended_signature = []
        for position in aligned_positions_extended_signature:
            extended_signature.append(aligned_domain[position])
        self.extended_signature = "".join(extended_signature)
