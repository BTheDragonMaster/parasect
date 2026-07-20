from typing import Optional
from dataclasses import dataclass
from enum import IntFlag
import logging

logger = logging.getLogger(__name__)

from Bio.SearchIO._model.hsp import HSP

@dataclass
class AOxHmm:
    ox_start: int = 215
    ox_end: int = 572
    ox_threshold: float = 0.8

    amp_upstream_start: int = 152
    amp_upstream_end: int = 214
    amp_downstream_start: int = 577
    amp_downstream_end: int = 751
    amp_threshold: float = 0.8

    def get_ox_length(self):
        return self.ox_end - self.ox_start

    def get_amp_length(self):
        return self.amp_upstream_end - self.amp_upstream_start + self.amp_downstream_end - self.amp_downstream_start

AOX_HMM = AOxHmm()

class DomainType(IntFlag):
    AMP_BINDING = 1
    A_OX = 2
    AMP_BINDING_C = 4

    N_TERMINAL = AMP_BINDING | A_OX

    def __str__(self):
        return self._to_string()

    def __repr__(self):
        return self._to_string()

    @classmethod
    def from_string(cls, string):
        from_string = {"AMP-binding": cls.AMP_BINDING,
                       "A-OX": cls.A_OX,
                       "AMP-binding_C": cls.AMP_BINDING_C}

        return from_string[string]

    def _to_string(self):
        to_string = {self.AMP_BINDING: "AMP-binding",
                     self.A_OX: "A-OX",
                     self.AMP_BINDING_C: "AMP-binding_C"}

        return to_string[self]


@dataclass
class HmmHit:
    """Class to store HMM hit"""
    protein_id: str
    domain_type: DomainType
    hsps: list[HSP]
    hmm_version: int

    def __repr__(self):
        return self._to_string()

    def __str__(self):
        return self._to_string()

    def __eq__(self, other):
        return self.__str__() == other.__str__()

    def _to_string(self):
        return f"{self.protein_id}|{self.domain_type}|{self.get_seq_start()}-{self.get_seq_end()}"

    def get_seq_start(self):
        return min(hsp.query_start for hsp in self.hsps)

    def get_seq_end(self):
        return max(hsp.query_end for hsp in self.hsps)

    def get_hmm_start(self):
        return min(hsp.hit_start for hsp in self.hsps)

    def get_hmm_end(self):
        return max(hsp.hit_end for hsp in self.hsps)


def _get_overlap_length(domain_start, domain_end, h_start, h_end):
    overlap = min([h_end, domain_end]) - max(domain_start, h_start)
    return max([overlap, 0])

def _resolve_n_terminal_hits(group: list[HmmHit]) -> list[HmmHit]:
    """Determine if the group of hits match best to an A-OX domain or an AMP-binding domain,
    only return those hits corresponding to the best-matching one. If only the OX domain is present,
    return an empty list.

    :param group: list of Hmm hits

    """

    has_amp = False
    has_ox = False
    ox_cover = 0
    amp_cover = 0

    for hit in group:
        if hit.domain_type == DomainType.A_OX:
            # Check if the A-OX HMM covers the OX-domain
            ox_cover += _get_overlap_length(hit.get_seq_start(),
                                            hit.get_seq_end(),
                                            AOX_HMM.ox_start,
                                            AOX_HMM.ox_end)

            # Check if the A-OX HMM covers the AMP-binding domain
            amp_cover += _get_overlap_length(hit.get_seq_start(),
                                             hit.get_seq_end(),
                                             AOX_HMM.amp_upstream_start,
                                             AOX_HMM.amp_upstream_end)

            amp_cover += _get_overlap_length(hit.get_seq_start(),
                                             hit.get_seq_end(),
                                             AOX_HMM.amp_downstream_start,
                                             AOX_HMM.amp_downstream_end)


    if ox_cover / AOX_HMM.get_ox_length() >= AOX_HMM.ox_threshold:
        logger.debug("A-OX domain found")
        has_ox = True

    if amp_cover / AOX_HMM.get_amp_length() >= AOX_HMM.amp_threshold:
        has_amp = True

    # Return the hits to the A-OX domain if it has both the OX and the AMP-binding domain
    if has_ox and has_amp:
        filtered_group = _filter_by_domain_type(group, DomainType.A_OX)
    # Return only the AMP-binding domain hits (possibly none, if it is an OX domain) otherwise
    else:
        filtered_group = _filter_by_domain_type(group, DomainType.AMP_BINDING)

    return filtered_group


def _group_hits(hits: list[HmmHit]) -> list[HmmHit]:
    if not hits:
        return []
    grouped_hits = []
    group = [hits[0]]

    for i, hit_1 in enumerate(hits):
        if i + 1 < len(hits):
            hit_2 = hits[i + 1]
            if hit_2.get_seq_start() - hit_1.get_seq_end() < 60:
                group.append(hit_2)
            else:
                grouped_hits.append(group[:])
                group = [hit_2]
        else:
            grouped_hits.append(group[:])
            group = []

    return grouped_hits

def _filter_by_domain_type(hits: list[HmmHit], domain_type: DomainType) -> list[HmmHit]:
    filtered_hits = []
    for hit in hits:
        if hit.domain_type == domain_type:
            filtered_hits.append(hit)

    return filtered_hits

def _merge_hits(hits: list[HmmHit]) -> Optional[HmmHit]:
    """
    Merge N-terminal AMP-binding hits

    :param hits: list of AMP-binding HMM hits
    """

    if hits:
        protein_id = hits[0].protein_id
        domain_type = hits[0].domain_type
        hmm_version = hits[0].hmm_version
        hsps = []

        for hit in hits:
            hsps.extend(hit.hsps)
            if hit.protein_id != protein_id:
                raise ValueError(f"Cannot merge hits from different sequences: {hit.protein_id}, {protein_id}")
            if hit.domain_type != domain_type:
                raise ValueError(f"Cannot merge different domain types: {hit.domain_type}, {domain_type}")
            if hit.hmm_version != hmm_version:
                raise ValueError(f"Cannot merge hits from different HMM versions: {hit.hmm_version}, {hmm_version}")

        merged_hit = HmmHit(protein_id, domain_type, hsps, hmm_version)
        return merged_hit
    else:
        raise ValueError("Need at least one hit for merging.")


def group_n_terminal_hits(hit_list: list[HmmHit]) -> list[HmmHit]:
    """
    Group and merge N-terminal AMP-binding hits within a single protein

    :param hit_list: list of AMP-binding HMM hits
    """
    n_terminal_hits = []
    c_terminal_hits = []
    seq_ids: set[str] = set()
    hmmer_versions: set[int] = set()

    for hit in hit_list:
        seq_ids.add(hit.protein_id)
        hmmer_versions.add(hit.hmm_version)

        if hit.domain_type & DomainType.N_TERMINAL:
            n_terminal_hits.append(hit)
        elif hit.domain_type == DomainType.AMP_BINDING_C:
            c_terminal_hits.append(hit)
        else:
            raise ValueError(f"Unsupported domain type: {hit.domain_type}")

    if len(seq_ids) > 1:
        raise ValueError("Cannot group hits from multiple sequences!")
    if len(hmmer_versions) > 1:
        raise ValueError("Cannot group hits from multiple hmmer versions!")

    n_terminal_hits.sort(key=lambda x: x.get_seq_start())
    c_terminal_hits.sort(key=lambda x: x.get_seq_start())

    grouped_n_terminal_hits = _group_hits(n_terminal_hits)

    merged_n_terminal = []

    for group in grouped_n_terminal_hits:

        resolved_group = _resolve_n_terminal_hits(group)
        if resolved_group:
            merged_hit = _merge_hits(resolved_group)
            merged_n_terminal.append(merged_hit)

    return merged_n_terminal + c_terminal_hits