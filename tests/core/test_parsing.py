# -*- coding: utf-8 -*-

"""Unit tests for the parts of parasect.core.parsing that the `paras` pipeline
touches: parse_fasta_file (used by featurisation.get_domains) and
parse_model_metadata_file (used by retrain_models.model_needs_retraining, on the
prepare_model() path).

Out of scope (PARASECT-only or model-training-only, not reached by `paras`):
parse_taxonomy_file, parse_raw_taxonomy, parse_list, parse_pcs, parse_esm_embedding,
parse_parasect_data, parse_smiles_mapping, parse_substrate_list,
data_from_substrate_names, bitvector_from_smiles, iterate_over_dir.
"""

import os
import unittest
from tempfile import TemporaryDirectory

from parasect.core.constants import MODEL_METADATA_FILE, get_path
from parasect.core.models import ModelType
from parasect.core.parsing import parse_fasta_file, parse_model_metadata_file

DOMAINS_FASTA_FILE = get_path("database_files/domains.fasta")


class TestParseFastaFile(unittest.TestCase):
    """Grounded in the real packaged domains.fasta database file."""

    @classmethod
    def setUpClass(cls):
        # Derive expected values by reading the raw file directly (one header + one
        # sequence line per record in this file), independently of parse_fasta_file
        # itself, rather than hand-typing a sequence from memory.
        with open(DOMAINS_FASTA_FILE) as fo:
            cls.first_header = fo.readline().strip()[1:]
            cls.first_sequence = fo.readline().strip()
        with open(DOMAINS_FASTA_FILE) as fo:
            cls.expected_record_count = sum(1 for line in fo if line.startswith(">"))

    def test_parses_every_record_in_the_real_database_fasta(self):
        fasta_dict = parse_fasta_file(DOMAINS_FASTA_FILE)
        self.assertEqual(len(fasta_dict), self.expected_record_count)
        self.assertEqual(fasta_dict[self.first_header], self.first_sequence)

    def test_missing_file_raises_filenotfounderror(self):
        with self.assertRaises(FileNotFoundError):
            parse_fasta_file("/no/such/file/anywhere.fasta")


class TestParseFastaFileStructuralEdgeCases(unittest.TestCase):
    """Small synthetic FASTA content for structural edge cases (sequence content
    itself is irrelevant here, so real domain data isn't needed for these)."""

    def _parse(self, content):
        with TemporaryDirectory() as tmp_dir:
            path = os.path.join(tmp_dir, "test.fasta")
            with open(path, "w") as fo:
                fo.write(content)
            return parse_fasta_file(path)

    def test_sequence_wrapped_across_multiple_lines_is_joined(self):
        result = self._parse(">h1\nAAA\nCCC\n>h2\nGGG\n")
        self.assertEqual(result, {"h1": "AAACCC", "h2": "GGG"})

    def test_empty_file_returns_empty_dict(self):
        self.assertEqual(self._parse(""), {})

    def test_trailing_blank_lines_do_not_corrupt_the_last_sequence(self):
        result = self._parse(">h1\nAAA\n\n\n")
        self.assertEqual(result, {"h1": "AAA"})

    def test_header_with_no_sequence_lines_is_silently_dropped_unless_last(self):
        # FLAG (possible bug): a header immediately followed by another header
        # (i.e. zero sequence lines) is dropped entirely instead of being kept
        # with an empty sequence -- *unless* it happens to be the final header in
        # the file, in which case it IS kept (with sequence ""). This asymmetry
        # means "h1" below vanishes silently rather than erroring or appearing
        # with an empty sequence. Flagging for Barbara to decide whether malformed
        # input like this should raise instead of losing a record silently.
        result = self._parse(">h1\n>h2\nAAA\n>h3\n")
        self.assertEqual(result, {"h2": "AAA", "h3": ""})
        self.assertNotIn("h1", result)


class TestParseModelMetadataFile(unittest.TestCase):
    def test_parses_the_real_packaged_metadata_file(self):
        model_to_version = parse_model_metadata_file(MODEL_METADATA_FILE)
        # The four models paras/parasect currently ship must all have an entry.
        for model_type in (ModelType.PARAS, ModelType.PARASECT, ModelType.PARAS_ALL_SUBSTRATES, ModelType.PARASECT_BACTERIAL):
            self.assertIn(model_type, model_to_version)
            self.assertIsInstance(model_to_version[model_type], str)

    def test_ignores_lines_not_starting_with_sklearn_version(self):
        with TemporaryDirectory() as tmp_dir:
            path = os.path.join(tmp_dir, "metadata.txt")
            with open(path, "w") as fo:
                fo.write("# a comment line\nsklearn_version_paras\t1.2.3\n")

            result = parse_model_metadata_file(path)

            self.assertEqual(result, {ModelType.PARAS: "1.2.3"})

    def test_empty_file_returns_empty_dict(self):
        with TemporaryDirectory() as tmp_dir:
            path = os.path.join(tmp_dir, "metadata.txt")
            with open(path, "w") as fo:
                fo.write("")

            self.assertEqual(parse_model_metadata_file(path), {})


if __name__ == "__main__":
    unittest.main()
