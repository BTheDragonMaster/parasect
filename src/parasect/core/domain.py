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
    new_positions = []
    position = 0

    # iterate over the aligned reference sequence
    for amino_acid_idx, amino_acid_id in enumerate(aligned_reference):
        if amino_acid_id != "-":
            if position in positions:
                new_positions.append(amino_acid_idx)

            position += 1

    return new_positions


def _get_reference_positions_hmm(
    query_sequence: str, reference_sequence: str, reference_positions: List[int]
) -> Optional[list[int]]:
    """Extract the given positions from a query alignment.

    :param query_sequence: The aligned query sequence.
    :type query_sequence: str
    :param reference_sequence: The aligned reference sequence.
    :type reference_sequence: str
    :param reference_positions: The positions of interest in the unaligned reference.
    :type reference_positions: List[int]
    :return: Positions of sequence elements at the adjusted reference,
        or None if the positions could not be extracted.
    :rtype: list[int
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
    return positions


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


def _get_gap_adjusted_positions(query: str, positions: list[int], offset: int) -> list[Optional[int]]:
    """Return sequence positions adjusted for gaps

    :param query: query sequence
    :type query: str
    :param positions: positions in the gapped query sequence
    :type positions: list[int]
    :param offset: query start position
    :type offset: int
    :return: list of gap-adjusted positions
    :rtype: list[int]
    """
    adjusted_positions: list[Optional[int]] = []
    position = offset
    for i, char in enumerate(query):
        if i in positions:
            if char != '-':
                adjusted_positions.append(position)
            else:
                adjusted_positions.append(None)
        if char != '-':
            position += 1

    return adjusted_positions


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
        self.domain_nr = 0
        self.start = domain_start
        self.end = domain_end

        self.sequence = ""
        self.protein_sequence = ""
        self.signature = ""
        self.extended_signature = ""
        self.signature_positions: list[int] = []
        self.extended_signature_positions: list[int] = []

    def domains_overlap(self, other: "AdenylationDomain", threshold: int = 50) -> bool:
        """
        Check if two domains overlap by at least a certain threshold of base pairs.

        :param other: Other adenylation domain
        :type other: AdenylationDomain
        :param threshold: The number of base pairs the two domains need to overlap by for a match
        :type threshold: int

        :return: bool, True if domains overlap, False otherwise
        """
        if other.start <= self.start <= other.end:
            if other.end - self.start >= threshold:
                return True
        if self.start <= other.start <= self.end:
            if self.end - other.start >= threshold:
                return True
        return False

    def set_domain_number(self, domain_nr: int) -> None:
        """Set the domain number.

        :param domain_nr: The domain number.
        :type domain_nr: int

        .. note:: This function modifies the domain number attribute.
        """
        self.domain_nr = domain_nr

    def set_protein_sequence(self, protein_sequence: str) -> None:
        """Set the protein sequence that contains the domain

        :param protein_sequence: The sequence of the domain.
        :type protein_sequence: str

        .. note:: This function modifies the sequence attribute.
        """
        self.protein_sequence = protein_sequence

    def set_sequence(self, sequence: str) -> None:
        """Set the sequence of the domain.

        :param sequence: The sequence of the domain.
        :type sequence: str

        .. note:: This function modifies the sequence attribute.
        """
        self.sequence = sequence

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
        query_offset = hit_n_terminal.query_start

        signature_location = _get_reference_positions_hmm(
            query_sequence=query,
            reference_sequence=profile,
            reference_positions=[p - offset for p in signature_positions],
        )

        if signature_location:
            signature = "".join([query[i] for i in signature_location])
            if all([char in valid for char in signature]):
                self.signature = signature
                self.signature_positions = _get_gap_adjusted_positions(query, signature_location, query_offset)

        lysine = None
        lysine_position = None
        query_c = None
        if hit_c_terminal:
            profile_c = hit_c_terminal.aln[1].seq
            query_c = hit_c_terminal.aln[0].seq
            offset_c = hit_c_terminal.hit_start

            lysine_position = _get_reference_positions_hmm(
                query_sequence=query_c,
                reference_sequence=profile_c,
                reference_positions=[p - offset_c for p in position_k],
            )

        if self.signature:
            if lysine_position and query_c:
                lysine = query_c[lysine_position[0]]
            if lysine and lysine in valid and lysine != "-":
                self.signature += lysine
            else:
                self.signature += "K"

        extended_signature_location = _get_reference_positions_hmm(
            query_sequence=query,
            reference_sequence=profile,
            reference_positions=[p - offset for p in extended_signature_positions],
        )
        if extended_signature_location:
            extended_signature = "".join([query[i] for i in extended_signature_location])
            if all([char in valid for char in extended_signature]):
                self.extended_signature = extended_signature
                self.extended_signature_positions = _get_gap_adjusted_positions(query, extended_signature_location,
                                                                                query_offset)

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
        self.signature_positions = _get_gap_adjusted_positions(aligned_domain, aligned_positions_signature, self.start)

        extended_signature = []
        for position in aligned_positions_extended_signature:
            extended_signature.append(aligned_domain[position])
        self.extended_signature = "".join(extended_signature)
        self.extended_signature_positions = _get_gap_adjusted_positions(aligned_domain,
                                                                        aligned_positions_extended_signature, self.start)
