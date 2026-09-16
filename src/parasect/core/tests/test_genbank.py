# -*- coding: utf-8 -*-

"""Unit tests for parasect.core.genbank (genbank_to_fasta, fetch_from_genbank)."""

import os
import unittest
from tempfile import TemporaryDirectory
from unittest.mock import patch
import subprocess

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from ncbi_acc_download.errors import DownloadError

from parasect.core.genbank import fetch_from_genbank, genbank_to_fasta


def write_genbank_record(path: str, cds_features: list) -> None:
    """Write a minimal, structurally real GenBank file via Biopython's own writer.

    :param cds_features: list of (start, end, qualifiers_dict) tuples, one per CDS feature.
    """
    record = SeqRecord(
        Seq("A" * 90), id="TEST_LOCUS", name="TEST_LOCUS", description="test record",
        annotations={"molecule_type": "DNA"},
    )
    for start, end, qualifiers in cds_features:
        record.features.append(SeqFeature(FeatureLocation(start, end), type="CDS", qualifiers=qualifiers))
    SeqIO.write(record, path, "genbank")


class TestGenbankToFasta(unittest.TestCase):
    def test_extracts_translations_keyed_by_protein_id_when_present(self):
        with TemporaryDirectory() as tmp_dir:
            in_path = os.path.join(tmp_dir, "in.gbk")
            out_path = os.path.join(tmp_dir, "out.fasta")
            write_genbank_record(in_path, [
                (0, 30, {"protein_id": ["PROT_1.1"], "translation": ["MEIKQACDEFGHIKLMNPQRSTVWYACDE"]}),
            ])

            genbank_to_fasta(in_path, out_path)

            with open(out_path) as fo:
                self.assertEqual(fo.read(), ">PROT_1.1\nMEIKQACDEFGHIKLMNPQRSTVWYACDE\n")

    def test_falls_back_to_locus_tag_when_protein_id_absent(self):
        with TemporaryDirectory() as tmp_dir:
            in_path = os.path.join(tmp_dir, "in.gbk")
            out_path = os.path.join(tmp_dir, "out.fasta")
            write_genbank_record(in_path, [
                (0, 30, {"locus_tag": ["locusA"], "translation": ["FGHIKLMNPQRSTVWYACDEFGHIKLMNP"]}),
            ])

            genbank_to_fasta(in_path, out_path)

            with open(out_path) as fo:
                self.assertEqual(fo.read(), ">locusA\nFGHIKLMNPQRSTVWYACDEFGHIKLMNP\n")

    def test_skips_cds_features_without_a_translation(self):
        with TemporaryDirectory() as tmp_dir:
            in_path = os.path.join(tmp_dir, "in.gbk")
            out_path = os.path.join(tmp_dir, "out.fasta")
            write_genbank_record(in_path, [
                (0, 30, {"gene": ["geneX"]}),  # no translation qualifier at all
            ])

            genbank_to_fasta(in_path, out_path)

            with open(out_path) as fo:
                self.assertEqual(fo.read(), "")

    def test_generates_gene_id_when_no_identifying_qualifier_present(self):
        with TemporaryDirectory() as tmp_dir:
            in_path = os.path.join(tmp_dir, "in.gbk")
            out_path = os.path.join(tmp_dir, "out.fasta")
            write_genbank_record(in_path, [
                (0, 30, {"translation": ["MEIKQACDEFGHIKLMNPQRSTVWYACDE"]}),
            ])

            genbank_to_fasta(in_path, out_path)

            with open(out_path) as fo:
                self.assertEqual(fo.read(), ">gene_0\nMEIKQACDEFGHIKLMNPQRSTVWYACDE\n")

    def test_missing_input_file_raises_filenotfounderror(self):
        with TemporaryDirectory() as tmp_dir:
            with self.assertRaises(FileNotFoundError):
                genbank_to_fasta(os.path.join(tmp_dir, "does_not_exist.gbk"), os.path.join(tmp_dir, "out.fasta"))


class TestFetchFromGenbank(unittest.TestCase):
    """fetch_from_genbank shells out to ncbi-acc-download, which needs live NCBI
    network access -- that part is not something a unit test should depend on, so
    subprocess.check_call is mocked here. An actual successful fetch against NCBI
    remains a manual/integration check.
    """

    @patch("parasect.core.genbank.subprocess.check_call")
    def test_builds_expected_command_for_multiple_accessions(self, mock_check_call):
        fetch_from_genbank(["ACC1", "ACC2"], "/tmp/out.fasta")

        mock_check_call.assert_called_once_with(
            ["ncbi-acc-download", "--format", "fasta", "--molecule", "protein", "--out", "/tmp/out.fasta", "ACC1", "ACC2"]
        )

    @patch("parasect.core.genbank.subprocess.check_call")
    def test_wraps_subprocess_failure_in_downloaderror(self, mock_check_call):
        mock_check_call.side_effect = subprocess.CalledProcessError(returncode=1, cmd="ncbi-acc-download")

        with self.assertRaises(DownloadError) as error:
            fetch_from_genbank(["BAD_ACCESSION"], "/tmp/out.fasta")

        self.assertEqual(str(error.exception), "Could not find one or more NCBI accessions")


if __name__ == "__main__":
    unittest.main()
