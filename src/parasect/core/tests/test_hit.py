import unittest
from unittest.mock import MagicMock

from Bio.SearchIO._model.hsp import HSP
from parasect.core.hit import _merge_hits, DomainType, HmmHit, group_n_terminal_hits


def make_hsp_mock(query_start: int,
                  query_end: int,
                  hit_start: int = 0,
                  hit_end: int = 0) -> MagicMock:
    """Mock for BioPython Bio.SearchIO.HSP class"""
    hsp = MagicMock(spec=HSP)
    hsp.query_start = query_start
    hsp.query_end = query_end
    hsp.hit_start = hit_start
    hsp.hit_end = hit_end
    return hsp

def make_hmm_hit(query_start: int,
                 query_end: int,
                 domain_type: DomainType = DomainType.AMP_BINDING,
                 hit_start: int = 0,
                 hit_end: int = 0,
                 protein_id: str = "mock_protein",
                 hmm_version: int = 3) -> HmmHit:
    hsp = make_hsp_mock(query_start, query_end, hit_start, hit_end)
    return HmmHit(protein_id=protein_id,
                  domain_type=domain_type,
                  hsps=[hsp],
                  hmm_version=hmm_version)

class TestHit(unittest.TestCase):
    hit_1 = make_hmm_hit(0, 165)
    hit_2 = make_hmm_hit(166, 300)
    hit_3 = make_hmm_hit(0, 165, domain_type=DomainType.AMP_BINDING_C)
    hit_4 = make_hmm_hit(167, 200, protein_id="mock_protein_2")
    hit_5 = make_hmm_hit(520, 600, domain_type=DomainType.AMP_BINDING_C)
    hit_6 = make_hmm_hit(560, 900)
    hit_7 = make_hmm_hit(359, 500)

    merged_hit_1 = make_hmm_hit(0, 300)
    merged_hit_2 = make_hmm_hit(0, 500)
    merged_hit_3 = make_hmm_hit(560, 900)

    def test_merge_hits(self):

        self.assertEqual(_merge_hits([self.hit_1, self.hit_2]), self.merged_hit_1)

        with self.assertRaises(ValueError):
            _merge_hits([self.hit_1, self.hit_3])

        with self.assertRaises(ValueError):
            _merge_hits([self.hit_1, self.hit_4])

        with self.assertRaises(ValueError):
            _merge_hits([])

    def test_group_n_terminal_hits(self):

        self.assertEqual(group_n_terminal_hits([self.hit_1, self.hit_2, self.hit_5, self.hit_6, self.hit_7]),
                         [self.merged_hit_2, self.merged_hit_3, self.hit_5])

        with self.assertRaises(ValueError):
            group_n_terminal_hits([self.hit_1, self.hit_2, self.hit_3, self.hit_4, self.hit_5, self.hit_6])


if __name__ == "__main__":
    unittest.main()
