# -*- coding: utf-8 -*-

"""Unit tests for parasect.core.tabular (Tabular, write_tabular).

Happy-path and column/row lookups are grounded in the real packaged
database_files/smiles.tsv file rather than a hand-rolled fixture.
"""

import os
import unittest
from tempfile import TemporaryDirectory

from parasect.core.constants import SMILES_FILE
from parasect.core.tabular import Tabular, write_tabular


class TestTabular(unittest.TestCase):
    """Tests against the real packaged smiles.tsv (substrate -> smiles)."""

    @classmethod
    def setUpClass(cls):
        cls.data = Tabular(path_in=SMILES_FILE, separator="\t")
        # Read the same row directly from disk so the expected value is never
        # hand-typed -- it comes from the real file, same as the class under test.
        with open(SMILES_FILE) as fo:
            header = fo.readline().strip().split("\t")
            first_data_line = fo.readline().strip().split("\t")
        cls.expected_columns = header
        cls.first_row_id = first_data_line[0]
        cls.first_row_smiles = first_data_line[1]

    def test_column_names_match_file_header(self):
        self.assertEqual(self.data.column_names, self.expected_columns)

    def test_get_row_value_matches_file_contents(self):
        self.assertEqual(self.data.get_row_value(self.first_row_id, "smiles"), self.first_row_smiles)

    def test_get_column_values_includes_every_row(self):
        smiles_column = self.data.get_column_values("smiles")
        self.assertEqual(len(smiles_column), len(self.data.rows))
        self.assertIn(self.first_row_smiles, smiles_column)

    def test_get_row_values_returns_full_row_in_column_order(self):
        self.assertEqual(self.data.get_row_values(self.first_row_id), [self.first_row_id, self.first_row_smiles])

    def test_unknown_column_raises_keyerror(self):
        with self.assertRaises(KeyError):
            self.data.get_column_values("not_a_real_column")

    def test_unknown_row_id_raises_keyerror(self):
        with self.assertRaises(KeyError):
            self.data.get_row_values("not_a_real_substrate_name")

    def test_missing_file_raises_filenotfounderror(self):
        with self.assertRaises(FileNotFoundError):
            Tabular(path_in="/no/such/file/anywhere.tsv")


class TestTabularMalformedInput(unittest.TestCase):
    """Structural edge cases for the tabular parser, using small synthetic files."""

    def _write(self, tmp_dir, name, content):
        path = os.path.join(tmp_dir, name)
        with open(path, "w") as fo:
            fo.write(content)
        return path

    def test_duplicate_row_id_raises_valueerror(self):
        with TemporaryDirectory() as tmp_dir:
            path = self._write(tmp_dir, "dup.tsv", "id\tvalue\nA\t1\nA\t2\n")
            with self.assertRaises(ValueError):
                Tabular(path_in=path)

    def test_row_with_wrong_number_of_columns_raises_valueerror(self):
        with TemporaryDirectory() as tmp_dir:
            path = self._write(tmp_dir, "ragged.tsv", "id\tvalue\nA\t1\nB\t1\t2\n")
            with self.assertRaises(ValueError):
                Tabular(path_in=path)

    def test_values_are_stripped_of_surrounding_whitespace(self):
        with TemporaryDirectory() as tmp_dir:
            path = self._write(tmp_dir, "whitespace.tsv", "id\tvalue\n A \t 1 \n")
            data = Tabular(path_in=path)
            self.assertEqual(data.get_row_value("A", "value"), "1")


class TestWriteTabular(unittest.TestCase):
    """Tests for the write_tabular helper."""

    def test_writes_header_and_rows_sorted_by_key(self):
        with TemporaryDirectory() as tmp_dir:
            out_file = os.path.join(tmp_dir, "out.tsv")
            write_tabular(
                dictionaries=[{"b_gene": "0.5", "a_gene": "0.1"}],
                header=["gene", "score"],
                out_file=out_file,
            )
            with open(out_file) as fo:
                content = fo.read()

            self.assertEqual(content, "gene\tscore\na_gene\t0.1\nb_gene\t0.5\n")

    def test_mismatched_header_and_dictionary_count_raises_assertionerror(self):
        with TemporaryDirectory() as tmp_dir:
            out_file = os.path.join(tmp_dir, "out.tsv")
            with self.assertRaises(AssertionError):
                write_tabular(dictionaries=[{"a": "1"}], header=["id", "col1", "col2"], out_file=out_file)

    def test_non_string_values_raise_typeerror(self):
        # FLAG: write_tabular builds each row with '\t'.join(row), which requires every
        # dictionary value to already be a str. Passing numeric values (a very plausible
        # real usage, e.g. writing out numeric scores) currently raises TypeError instead
        # of being coerced. Documenting the current behavior here rather than working
        # around it -- flagging for Barbara to decide whether values should be str()'d.
        with TemporaryDirectory() as tmp_dir:
            out_file = os.path.join(tmp_dir, "out.tsv")
            with self.assertRaises(TypeError):
                write_tabular(dictionaries=[{"a_gene": 0.1}], header=["gene", "score"], out_file=out_file)


if __name__ == "__main__":
    unittest.main()
