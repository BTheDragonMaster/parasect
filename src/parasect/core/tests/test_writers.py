# -*- coding: utf-8 -*-

"""Unit tests for parasect.core.writers."""

import os
import unittest
from tempfile import TemporaryDirectory

from parasect.core.domain import AdenylationDomain
from parasect.core.hit import DomainType
from parasect.core.models import ModelType
from parasect.core.parasect_result import Result
from parasect.core.parsing import parse_model_metadata_file
from parasect.core.writers import write_fasta_file, write_list, write_model_metadata_file, write_results


def make_domain(protein_name: str, domain_nr: int, start: int, end: int,
                sequence: str, signature: str, extended_signature: str) -> AdenylationDomain:
    domain = AdenylationDomain(protein_name, DomainType.AMP_BINDING, start, end)
    domain.set_domain_number(domain_nr)
    domain.set_sequence(sequence)
    domain.signature = signature
    domain.extended_signature = extended_signature
    return domain


class TestWriteFastaFile(unittest.TestCase):
    def test_writes_headers_in_sorted_order(self):
        with TemporaryDirectory() as tmp_dir:
            out_path = os.path.join(tmp_dir, "out.fasta")
            write_fasta_file({"zebra": "AAA", "alpha": "CCC"}, out_path)
            with open(out_path) as fo:
                content = fo.read()

            self.assertEqual(content, ">alpha\nCCC\n>zebra\nAAA\n")

    def test_empty_dict_writes_empty_file(self):
        with TemporaryDirectory() as tmp_dir:
            out_path = os.path.join(tmp_dir, "empty.fasta")
            write_fasta_file({}, out_path)
            with open(out_path) as fo:
                content = fo.read()

            self.assertEqual(content, "")


class TestWriteList(unittest.TestCase):
    def test_sorts_by_default(self):
        with TemporaryDirectory() as tmp_dir:
            out_path = os.path.join(tmp_dir, "list.txt")
            write_list(["banana", "apple", "cherry"], out_path)
            with open(out_path) as fo:
                lines = fo.read().splitlines()

            self.assertEqual(lines, ["apple", "banana", "cherry"])

    def test_preserves_order_when_sort_is_false(self):
        with TemporaryDirectory() as tmp_dir:
            out_path = os.path.join(tmp_dir, "list.txt")
            write_list(["banana", "apple", "cherry"], out_path, sort=False)
            with open(out_path) as fo:
                lines = fo.read().splitlines()

            self.assertEqual(lines, ["banana", "apple", "cherry"])

    def test_empty_list_writes_empty_file(self):
        with TemporaryDirectory() as tmp_dir:
            out_path = os.path.join(tmp_dir, "list.txt")
            write_list([], out_path)
            with open(out_path) as fo:
                content = fo.read()

            self.assertEqual(content, "")


class TestWriteModelMetadataFile(unittest.TestCase):
    def test_round_trips_through_the_real_parser(self):
        # Rather than hand-writing the expected file text, verify the round trip
        # through parasect's own parser -- if either side's format assumptions
        # drift, this breaks instead of silently reading back the wrong thing.
        model_to_version = {ModelType.PARAS: "1.2.0", ModelType.PARASECT: "1.3.4"}
        with TemporaryDirectory() as tmp_dir:
            out_path = os.path.join(tmp_dir, "metadata.txt")
            write_model_metadata_file(model_to_version, out_path)

            parsed_back = parse_model_metadata_file(out_path)

        self.assertEqual(parsed_back, model_to_version)


class TestWriteResults(unittest.TestCase):
    def _make_results(self):
        domain_1 = make_domain("dptA", 1, 10, 175, "MEIKQ", "ACDEFGHIKL", "ACDEFGHIKLMNPQRSTVWY")
        domain_2 = make_domain("dptB", 1, 5, 160, "QKIEM", "LKIHGFEDCA", "YWVTSRQPNMLKIHGFEDCA")

        result_1 = Result(domain_1, [0.1, 0.7, 0.2], ["alanine", "valine", "glycine"],
                          ["CC(N)C(=O)O", "CC(C)C(N)C(=O)O", "NCC(=O)O"])
        result_2 = Result(domain_2, [0.9, 0.05, 0.05], ["glycine", "alanine", "valine"],
                          ["NCC(=O)O", "CC(N)C(=O)O", "CC(C)C(N)C(=O)O"])
        return [result_1, result_2]

    def test_raises_if_more_predictions_requested_than_available(self):
        results = self._make_results()
        with TemporaryDirectory() as tmp_dir:
            with self.assertRaises(ValueError):
                write_results(results, tmp_dir, number_predictions=10, model_type=ModelType.PARAS)

    def test_main_results_file_reports_top_n_sorted_predictions(self):
        results = self._make_results()
        with TemporaryDirectory() as tmp_dir:
            write_results(results, tmp_dir, number_predictions=2, model_type=ModelType.PARAS, job_name="myjob")

            result_file = os.path.join(tmp_dir, "myjob_paras_results.txt")
            self.assertTrue(os.path.exists(result_file))
            with open(result_file) as fo:
                lines = fo.read().splitlines()

            self.assertEqual(lines[0], "domain_id\tprediction_1\tconfidence_prediction_1\tprediction_2\tconfidence_prediction_2")
            self.assertEqual(lines[1], "dptA|domain_1|10-175\tvaline\t0.7\tglycine\t0.2")
            self.assertEqual(lines[2], "dptB|domain_1|5-160\tglycine\t0.9\talanine\t0.05")

    def test_write_results_mutates_results_in_place_via_sort(self):
        # FLAG: write_results() calls result.sort() on every Result it's given,
        # mutating the caller's objects as a side effect rather than working on
        # copies. Documenting this rather than changing it.
        results = self._make_results()
        self.assertEqual(results[0].predictions, [0.1, 0.7, 0.2])

        with TemporaryDirectory() as tmp_dir:
            write_results(results, tmp_dir, number_predictions=1, model_type=ModelType.PARAS)

        self.assertEqual(results[0].predictions, [0.7, 0.2, 0.1])

    def test_optional_fasta_outputs_are_only_written_when_requested(self):
        results = self._make_results()
        with TemporaryDirectory() as tmp_dir:
            write_results(results, tmp_dir, number_predictions=1, model_type=ModelType.PARAS, job_name="myjob")

            for suffix in ("signatures", "extended_signatures", "sequences"):
                self.assertFalse(os.path.exists(os.path.join(tmp_dir, f"myjob_{suffix}.fasta")))

    def test_optional_fasta_outputs_contain_correct_content_when_requested(self):
        results = self._make_results()
        with TemporaryDirectory() as tmp_dir:
            write_results(
                results, tmp_dir, number_predictions=1, model_type=ModelType.PARAS, job_name="myjob",
                save_signatures=True, save_extended_signatures=True, save_domain_sequences=True,
            )

            with open(os.path.join(tmp_dir, "myjob_signatures.fasta")) as fo:
                self.assertEqual(fo.read(), ">dptA|domain_1|10-175\nACDEFGHIKL\n>dptB|domain_1|5-160\nLKIHGFEDCA\n")

            with open(os.path.join(tmp_dir, "myjob_sequences.fasta")) as fo:
                self.assertEqual(fo.read(), ">dptA|domain_1|10-175\nMEIKQ\n>dptB|domain_1|5-160\nQKIEM\n")


if __name__ == "__main__":
    unittest.main()
