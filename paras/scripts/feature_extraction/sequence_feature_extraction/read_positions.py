from typing import List
import os

import paras.data.sequence_data.residue_positions

APOSITION_FILE = os.path.join(os.path.dirname(paras.data.sequence_data.residue_positions.__file__), 'stachelhaus.txt')
APOSITION_FILE_34 = os.path.join(os.path.dirname(paras.data.sequence_data.residue_positions.__file__), 'active_site.txt')
APOSITION_FILE_HMM3 = os.path.join(os.path.dirname(paras.data.sequence_data.residue_positions.__file__), 'stachelhaus_hmm3.txt')
APOSITION_FILE_34_HMM3 = os.path.join(os.path.dirname(paras.data.sequence_data.residue_positions.__file__), 'active_site_hmm3.txt')
APOSITION_FILE_HMM2 = os.path.join(os.path.dirname(paras.data.sequence_data.residue_positions.__file__), 'stachelhaus_hmm2.txt')
APOSITION_FILE_34_HMM2 = os.path.join(os.path.dirname(paras.data.sequence_data.residue_positions.__file__), 'active_site_hmm2.txt')


START_POSITION = 66


def read_positions(position_file: str, start_position: int) -> List[int]:
    """
    Return positions from a tab-separated file. Positions are relative to start_position.

    Input:
    position_file: str, the path to tab-separated file containing the positions
    start_position: int, a relative start position to adjust all positions by

    Output:
    positions: list of [int, ->], one for each position found in the file
    """

    with open(position_file, 'r') as position_data:
        text = position_data.read().strip()
        positions = []
        for i in text.split("\t"):
            positions.append(int(i) - start_position)
    return positions


def get_reference_positions(positions: List[int], aligned_reference: str) -> List[int]:
    """
    Adjusts a list of positions to account for gaps in the reference sequence

    Input:
    positions: list of [int, ->], with integers representing positions of interest in
        the reference sequence
    aligned_reference: the (aligned) reference sequence

    Output:
    pos_list: list of [int, ->], a new list of positions, each >= the original position
    """
    pos_list = []
    position = 0
    for i, aa in enumerate(aligned_reference):
        if aa != "-":
            if position in positions:
                pos_list.append(i)
            position += 1
    return pos_list


# From antiSMASH 7.0
def get_reference_positions_hmm(query, reference, ref_positions):
    """ Extracts the given positions from a query alignment. The positions are
        adjusted to account for any gaps in the reference sequence.

        Arguments:
            query: the aligned query
            reference: the aligned reference
            ref_positions: the positions of interest in the unaligned reference

        Returns:
            a string containing the sequence elements at the adjusted reference
            positions or None if the reference is too short for some reason
    """
    # adjust position of interest to account for gaps in the ref sequence alignment
    positions = []
    position_skipping_gaps = 0
    for i, amino in enumerate(reference):
        if amino in "-.":
            continue
        if position_skipping_gaps in ref_positions:
            positions.append(i)
        position_skipping_gaps += 1
    if len(positions) != len(ref_positions):
        return None
    # extract positions from query sequence
    return "".join([query[i] for i in positions])


POSITIONS_SIGNATURE = read_positions(APOSITION_FILE, START_POSITION)
POSITIONS_EXTENDED_SIGNATURE = read_positions(APOSITION_FILE_34, START_POSITION)

HMM2_POSITIONS_SIGNATURE = read_positions(APOSITION_FILE_HMM2, 0)
HMM2_POSITIONS_EXTENDED_SIGNATURE = read_positions(APOSITION_FILE_34_HMM2, 0)

HMM3_POSITIONS_SIGNATURE = read_positions(APOSITION_FILE_HMM3, 0)
HMM3_POSITIONS_EXTENDED_SIGNATURE = read_positions(APOSITION_FILE_34_HMM3, 0)

# Hmmer2:
HMM2_POSITION_K = [36]

# Hmmer3:
HMM3_POSITION_K = [75]
