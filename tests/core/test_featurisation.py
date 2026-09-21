import os
import unittest
from math import isclose
from tempfile import TemporaryDirectory
from unittest.mock import patch

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from parasect.core.featurisation import get_domain_features, get_domains


def genbank_text(qualifiers: dict) -> str:
    """Return a minimal GenBank file with a single CDS feature, as text."""
    record = SeqRecord(Seq("A" * 30), id="TEST_LOCUS", name="TEST_LOCUS", annotations={"molecule_type": "DNA"})
    record.features.append(SeqFeature(FeatureLocation(0, 30), type="CDS", qualifiers=qualifiers))
    return record.format("genbank")


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


class TestGetDomainsWithoutSequences(unittest.TestCase):
    """Input yielding no protein sequences must fail before HMMER, naming the likely cause."""

    def _assert_no_sequences_error(self, content: str, suffix: str, file_type: str, expected: str) -> None:
        with TemporaryDirectory() as tmp_dir:
            path_in = os.path.join(tmp_dir, f"input.{suffix}")
            with open(path_in, "w") as fo:
                fo.write(content)

            with patch("parasect.core.featurisation.run_hmmpfam2") as run_hmmpfam2:
                with self.assertRaises(ValueError) as ctx:
                    get_domains(path_in, tmp_dir, "hmm", file_type)

            run_hmmpfam2.assert_not_called()
            self.assertIn(expected, str(ctx.exception))

    def test_genbank_parsed_as_fasta_names_mismatch(self):
        gbk = genbank_text({"protein_id": ["P1"], "translation": ["MEIKQ"]})
        self._assert_no_sequences_error(gbk, "gbk", "fasta", "Set the input type to GBK")

    def test_fasta_parsed_as_genbank_names_mismatch(self):
        self._assert_no_sequences_error(">seq_1\nMEIKQ\n", "fasta", "gbk", "Set the input type to FASTA")

    def test_genbank_without_translations(self):
        gbk = genbank_text({"protein_id": ["P1"]})
        self._assert_no_sequences_error(gbk, "gbk", "gbk", "no CDS features with a translation")

    def test_empty_fasta(self):
        self._assert_no_sequences_error("", "fasta", "fasta", "no protein sequences found in FASTA input")


if __name__ == "__main__":
    unittest.main()
