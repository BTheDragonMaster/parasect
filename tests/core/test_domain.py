import unittest
from parasect.core.domain import _merge_signatures


class TestSignatureMerging(unittest.TestCase):
    def test_merge_signatures(self):
        signatures = [['-', 'C', '-', 'E', 'F'],
                      ['A', 'C', 'D', '-', '-']]
        consensus_signature = "ACDEF"
        consensus_positions = [1, 2, 3, 4, 5]

        positions = [[None, 2, None, 4, 5],
                     [1, 2, 3, None, None]]

        signature, signature_positions = _merge_signatures(signatures, positions)

        self.assertEqual(signature, consensus_signature)
        self.assertEqual(signature_positions, consensus_positions)
