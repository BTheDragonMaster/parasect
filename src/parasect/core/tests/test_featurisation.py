import unittest
from math import isclose
from parasect.core.featurisation import get_domain_features, merge_hits, group_n_terminal_hits


class TestFeaturisation(unittest.TestCase):
    def test_get_domain_features(self):
        sequence_1 = 'AAA'
        vector = [0.07, -1.73, 0.09, 0, 8.1, -0.06, 0.00, 90.0, 1.42, 0.83, 0.66, 6.00, 0.06, -0.25, 0.25,
                  0.07, -1.73, 0.09, 0, 8.1, -0.06, 0.00, 90.0, 1.42, 0.83, 0.66, 6.00, 0.06, -0.25, 0.25,
                  0.07, -1.73, 0.09, 0, 8.1, -0.06, 0.00, 90.0, 1.42, 0.83, 0.66, 6.00, 0.06, -0.25, 0.25]
        self.assertNearlyEqual(get_domain_features(sequence_1), vector)

    def test_merge_hits(self):
        hit_1 = ("AMP-binding", 0, 165, "seq_1|AMP-binding|0-165")
        hit_2 = ("AMP-binding", 166, 300, "seq_1|AMP-binding|166-300")
        hit_3 = ("AMP-binding_C", 0, 165, "seq_1|AMP-binding_C|0-165")
        hit_4 = ("AMP-binding", 167, 200, "seq_2|AMP-binding|167-200")

        self.assertEqual(merge_hits([hit_1, hit_2]), ("AMP-binding", 0, 300, "seq_1|AMP-binding|0-300"))

        with self.assertRaises(ValueError):
            merge_hits([hit_1, hit_3])

        with self.assertRaises(ValueError):
            merge_hits([hit_1, hit_4])

        with self.assertRaises(ValueError):
            merge_hits([])

    def test_group_n_terminal_hits(self):
        hit_1 = ("AMP-binding", 0, 165, "seq_1|AMP-binding|0-165")
        hit_2 = ("AMP-binding", 166, 300, "seq_1|AMP-binding|166-300")
        hit_3 = ("AMP-binding_C", 520, 600, "seq_1|AMP-binding_C|520-600")
        hit_4 = ("AMP-binding", 167, 200, "seq_2|AMP-binding|167-200")
        hit_5 = ("AMP-binding", 560, 900, "seq_1|AMP-binding|560-900")
        hit_6 = ("AMP-binding", 359, 500, "seq_1|AMP-binding|500-900")

        merged_hit_1 = ("AMP-binding", 0, 500, "seq_1|AMP-binding|0-500")
        merged_hit_2 = ("AMP-binding", 560, 900, "seq_1|AMP-binding|560-900")

        self.assertEqual(group_n_terminal_hits([hit_1, hit_2, hit_3, hit_5, hit_6]),
                         [merged_hit_1, merged_hit_2, hit_3])

        with self.assertRaises(ValueError):
            group_n_terminal_hits([hit_1, hit_2, hit_3, hit_4, hit_5, hit_6])

    def assertNearlyEqual(self, list_1, list_2):
        for i, element_1 in enumerate(list_1):
            element_2 = list_2[i]
            if not isclose(element_1, element_2, rel_tol=0.0001):
                self.fail(f"Lists are not equal: {list_1}, {list_2}. \n First mismatching element: {i} ([{element_1}], [{element_2}])")


if __name__ == "__main__":
    unittest.main()
