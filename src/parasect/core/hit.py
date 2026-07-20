from typing import Optional
from dataclasses import dataclass
from enum import Enum

from Bio.SearchIO._model.hsp import HSP

from parasect.core.constants import (OX_START_POSITION, OX_END_POSITION, AMP_UP_START_POSITION, AMP_UP_END_POSITION,
                                     AMP_DOWN_START_POSITION, AMP_DOWN_END_POSITION, OX_THRESHOLD,
                                     AMP_THRESHOLD, AMP_LENGTH, OX_LENGTH)


class HitType(Enum):
    AMP_BINDING = 1
    A_OX = 2
    AMP_BINDING_C = 3

    @classmethod
    def from_string(cls, string):
        from_string = {"AMP-binding": cls.AMP_BINDING,
                       "A-OX": cls.A_OX,
                       "AMP-binding_C": cls.AMP_BINDING_C}

        return from_string[string]

@dataclass
class HmmHit:
    """Class to store HMM hit"""
    id: str
    hsp: HSP
    hit_type: HitType
    hmm_version: int

    def get_seq_start(self):
        return self.hsp.query_start

    def get_seq_end(self):
        return self.hsp.query_end

    def get_hmm_start(self):
        return self.hsp.hit_start

    def get_hmm_end(self):
        return self.hsp.hit_end


def _get_overlap_length(domain_start, domain_end, h_start, h_end):
    overlap = min([h_end, domain_end]) - max(domain_start, h_start)
    return max([overlap, 0])

def _get_domain_type(group, hit_lookup):

    has_amp = False
    has_ox = False
    ox_cover = 0
    amp_cover = 0

    for ox_hit in group:
        full_hit = hit_lookup[ox_hit[3]]
        print(full_hit.hit_start, full_hit.hit_end, full_hit.id)
        ox_cover += _get_overlap_length(full_hit.hit_start,
                                        full_hit.hit_end,
                                        OX_START_POSITION,
                                        OX_END_POSITION)

        amp_cover += _get_overlap_length(full_hit.hit_start,
                                         full_hit.hit_end,
                                         AMP_UP_START_POSITION,
                                         AMP_UP_END_POSITION)

        amp_cover += _get_overlap_length(full_hit.hit_start,
                                         full_hit.hit_end,
                                         AMP_DOWN_START_POSITION,
                                         AMP_DOWN_END_POSITION)

    if ox_cover / OX_LENGTH >= OX_THRESHOLD:
        has_ox = True
    if amp_cover / AMP_LENGTH >= AMP_THRESHOLD:
        has_amp = True

    if has_ox and has_amp:
        return "A-OX"
    elif has_amp:
        return "AMP-binding"
    else:
        return None


def _group_hits(hits):
    if not hits:
        return []
    grouped_hits = []
    group = [hits[0]]

    for i, hit_1 in enumerate(hits):
        if i + 1 < len(hits):
            hit_2 = hits[i + 1]
            if hit_2[1] - hit_1[2] < 60:
                group.append(hit_2)
            else:
                grouped_hits.append(group[:])
                group = [hit_2]
        else:
            grouped_hits.append(group[:])
            group = []

    return grouped_hits

def _filter_by_domain_type(hits: list[tuple[str, int, int, str]], domain_type: str) -> list[tuple[str, int, int, str]]:
    if domain_type not in ["AMP-binding", "A-OX"]:
        raise ValueError(f"Unknown domain type: {domain_type}")
    filtered_hits = []
    for n_hit in hits:
        if n_hit[0] == domain_type:
            filtered_hits.append(n_hit)

    return filtered_hits

def merge_hits(hits: list[tuple[str, int, int, str]], domain_type: str) -> Optional[tuple[str, int, int, str]]:
    """
    Merge N-terminal AMP-binding hits

    :param hits: list of AMP-binding HMM hits
    :type hits: list[tuple[str, int, int, str]]
    :param domain_type: Type of N-terminal AMP-binding domain (AMP-binding or A-OX)
    """

    hits = _filter_by_domain_type(hits, domain_type)

    if hits:
        seq_id, hit_id, _ = hits[0][3].split('|')
        for hit in hits:
            seq_id_2, hit_id_2, _ = hit[3].split('|')
            if seq_id_2 != seq_id:
                raise ValueError(f"Cannot merge hits from different sequences! {seq_id}, {seq_id_2}")
            if hit_id_2 != hit_id:
                raise ValueError(f"Cannot merge different hit types! {hit_id}, {hit_id_2}")

        hit_start = min([hit[1] for hit in hits])
        hit_end = max([hit[2] for hit in hits])
        hit_key = f"{seq_id}|{hit_id}|{hit_start}-{hit_end}"
        merged_hit = (hit_id, hit_start, hit_end, hit_key)
        return merged_hit
    else:
        return None


def group_n_terminal_hits(hit_list: list[tuple[str, int, int, str]],
                          id_to_hit: dict[str, HSP]) -> tuple[list[tuple[str, int, int, str]], dict[str, list[str]]]:
    """
    Group and merge N-terminal AMP-binding hits within a single protein

    :param hit_list: list of AMP-binding HMM hits
    :type hit_list: list[tuple[str, int, int, str]]
    :param id_to_hit: dictionary of protein ID to HSPs within that protein
    """
    n_terminal_hits = []
    c_terminal_hits = []
    seq_ids = set()

    for hit in hit_list:

        hit_id, hit_start, hit_end, hit_key = hit
        seq_id = hit_key.split('|')[0]
        seq_ids.add(seq_id)

        if hit_id in ["AMP-binding", "A-OX"]:
            n_terminal_hits.append(hit)
        elif hit_id == "AMP-binding_C":
            c_terminal_hits.append(hit)

    if len(seq_ids) > 1:
        raise ValueError("Cannot group hits from multiple sequences!")

    n_terminal_hits.sort(key=lambda x: x[1])
    c_terminal_hits.sort(key=lambda x: x[1])

    grouped_n_terminal_hits = _group_hits(n_terminal_hits)

    merged_n_terminal = []
    merged_to_original = {}

    for n_group in grouped_n_terminal_hits:
        domain_type = _get_domain_type(n_group, id_to_hit)
        merged_hit = merge_hits(n_group, domain_type)
        merged_n_terminal.append(merged_hit)
        merged_to_original[merged_hit[3]] = []
        for hit in n_group:
            merged_to_original[merged_hit[3]].append(hit[3])

    return merged_n_terminal + c_terminal_hits, merged_to_original