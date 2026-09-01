import unittest
from math import isclose
from parasect.core.featurisation import get_domain_features, group_n_terminal_hits


class TestFeaturisation(unittest.TestCase):
    def test_get_domain_features(self):
        sequence_1 = 'AAA'
        vector = [0.07, -1.73, 0.09, 0, 8.1, -0.06, 0.00, 90.0, 1.42, 0.83, 0.66, 6.00, 0.06, -0.25, 0.25,
                  0.07, -1.73, 0.09, 0, 8.1, -0.06, 0.00, 90.0, 1.42, 0.83, 0.66, 6.00, 0.06, -0.25, 0.25,
                  0.07, -1.73, 0.09, 0, 8.1, -0.06, 0.00, 90.0, 1.42, 0.83, 0.66, 6.00, 0.06, -0.25, 0.25]
        self.assertNearlyEqual(get_domain_features(sequence_1), vector)


    def assertNearlyEqual(self, list_1, list_2):
        for i, element_1 in enumerate(list_1):
            element_2 = list_2[i]
            if not isclose(element_1, element_2, rel_tol=0.0001):
                self.fail(f"Lists are not equal: {list_1}, {list_2}. \n First mismatching element: {i} ([{element_1}], [{element_2}])")


if __name__ == "__main__":
    unittest.main()
