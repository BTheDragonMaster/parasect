# -*- coding: utf-8 -*-

"""Unit tests for parasect.core.parasect_result.Result."""

import unittest

from parasect.core.domain import AdenylationDomain
from parasect.core.hit import DomainType
from parasect.core.parasect_result import Result


def make_domain(protein_name: str = "dptA", domain_nr: int = 1, start: int = 10, end: int = 175,
                sequence: str = "MEIKQ", signature: str = "ACDEFGHIKL",
                extended_signature: str = "ACDEFGHIKLMNPQRSTVWY") -> AdenylationDomain:
    """Build a real (not mocked) AdenylationDomain for use as test fixture data."""
    domain = AdenylationDomain(protein_name, DomainType.AMP_BINDING, start, end)
    domain.set_domain_number(domain_nr)
    domain.set_sequence(sequence)
    domain.signature = signature
    domain.extended_signature = extended_signature
    return domain


class TestResult(unittest.TestCase):
    def test_sort_orders_predictions_labels_and_smiles_together_descending(self):
        result = Result(
            domain=make_domain(),
            predictions=[0.1, 0.7, 0.2],
            prediction_labels=["alanine", "valine", "glycine"],
            prediction_smiles=["CC(N)C(=O)O", "CC(C)C(N)C(=O)O", "NCC(=O)O"],
        )
        result.sort()

        self.assertEqual(result.predictions, [0.7, 0.2, 0.1])
        self.assertEqual(result.prediction_labels, ["valine", "glycine", "alanine"])
        # smiles are private but exposed through to_json; check they moved with their label
        json_predictions = result.to_json()["predictions"]
        self.assertEqual(
            [p["substrate_smiles"] for p in json_predictions],
            ["CC(C)C(N)C(=O)O", "NCC(=O)O", "CC(N)C(=O)O"],
        )

    def test_sort_is_stable_for_tied_predictions(self):
        # Python's list.sort() is stable: equal-probability entries should keep
        # their original relative order rather than being reshuffled.
        result = Result(
            domain=make_domain(),
            predictions=[0.5, 0.9, 0.5],
            prediction_labels=["first_tied", "top", "second_tied"],
            prediction_smiles=["C", "CC", "CCC"],
        )
        result.sort()

        self.assertEqual(result.predictions, [0.9, 0.5, 0.5])
        self.assertEqual(result.prediction_labels, ["top", "first_tied", "second_tied"])

    def test_get_domain_header_uses_default_separators(self):
        domain = make_domain(protein_name="dptA", domain_nr=2, start=10, end=175)
        result = Result(domain, [1.0], ["alanine"], ["CC(N)C(=O)O"])

        self.assertEqual(result.get_domain_header(), "dptA|domain_2|10-175")

    def test_get_domain_header_uses_custom_separators(self):
        domain = make_domain(protein_name="dptA", domain_nr=2, start=10, end=175)
        result = Result(domain, [1.0], ["alanine"], ["CC(N)C(=O)O"])

        self.assertEqual(result.get_domain_header(separator_1=":", separator_2="#", separator_3="~"), "dptA:domain#2:10~175")

    def test_to_json_reports_domain_fields_and_zipped_predictions(self):
        domain = make_domain(
            protein_name="dptA", domain_nr=1, start=10, end=175,
            sequence="MEIKQ", signature="ACDEFGHIKL", extended_signature="ACDEFGHIKLMNPQRSTVWY",
        )
        result = Result(
            domain,
            predictions=[0.3, 0.6],
            prediction_labels=["alanine", "valine"],
            prediction_smiles=["CC(N)C(=O)O", "CC(C)C(N)C(=O)O"],
        )

        data = result.to_json()

        self.assertEqual(data["domain_name"], "dptA")
        self.assertEqual(data["domain_nr"], 1)
        self.assertEqual(data["domain_start"], 10)
        self.assertEqual(data["domain_end"], 175)
        self.assertEqual(data["domain_sequence"], "MEIKQ")
        self.assertEqual(data["domain_signature"], "ACDEFGHIKL")
        self.assertEqual(data["domain_extended_signature"], "ACDEFGHIKLMNPQRSTVWY")
        self.assertEqual(
            data["predictions"],
            [
                {"substrate_name": "alanine", "substrate_smiles": "CC(N)C(=O)O", "probability": 0.3},
                {"substrate_name": "valine", "substrate_smiles": "CC(C)C(N)C(=O)O", "probability": 0.6},
            ],
        )

    def test_predictions_and_prediction_labels_properties_reflect_state_after_sort(self):
        result = Result(make_domain(), [0.2, 0.8], ["a", "b"], ["C", "CC"])
        self.assertEqual(result.predictions, [0.2, 0.8])
        self.assertEqual(result.prediction_labels, ["a", "b"])

        result.sort()

        self.assertEqual(result.predictions, [0.8, 0.2])
        self.assertEqual(result.prediction_labels, ["b", "a"])


if __name__ == "__main__":
    unittest.main()
