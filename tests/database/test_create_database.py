"""Fresh database creation includes all input data and ordered substrate links."""

import sqlite3
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from sqlalchemy import create_engine, select
from sqlalchemy.orm import Session

from parasect.database.build_database import AdenylationDomain
from parasect.database.create_database import create_database


class TestCreateDatabase(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.data_dir = self.root / "inputs"
        self.data_dir.mkdir()
        self.output = self.root / "new.db"
        (self.data_dir / "parasect_dataset.txt").write_text(
            "domain_id\tsequence\tspecificity\nTest.A1\tAAAA\tvaline|alanine|leucine\n",
        )
        (self.data_dir / "signatures.fasta").write_text(">Test.A1\n" + "A" * 10 + "\n")
        (self.data_dir / "extended_signatures.fasta").write_text(">Test.A1\n" + "A" * 34 + "\n")
        (self.data_dir / "proteins.fasta").write_text(">Test\nMAAAAK\n")
        (self.data_dir / "smiles.tsv").write_text(
            "substrate\tsmiles\n"
            "alanine\tC[C@H](N)C(=O)O\n"
            "leucine\tCC(C)C[C@H](N)C(=O)O\n"
            "valine\tCC(C)[C@H](N)C(=O)O\n",
        )
        (self.data_dir / "taxonomy.txt").write_text(
            "protein_id\tdomain\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tstrain\n"
            "Test\tBacteria\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tTest\tTest species\tUnknown\n",
        )

    def test_creates_populated_database_with_persistent_substrate_order(self):
        create_database(self.output, self.data_dir)
        engine = create_engine(f"sqlite:///{self.output}")
        self.addCleanup(engine.dispose)
        with Session(engine) as session:
            domain = session.scalars(select(AdenylationDomain)).one()
            self.assertEqual([s.name for s in domain.substrates], ["valine", "alanine", "leucine"])
            self.assertEqual([a.position for a in domain.substrate_associations], [0, 1, 2])
            self.assertEqual(domain.get_name(), "Test.A1")
            self.assertEqual(domain.proteins[0].protein.get_name(), "Test")
            self.assertEqual(domain.proteins[0].protein.taxonomy.species, "Test species")
        connection = sqlite3.connect(self.output)
        self.addCleanup(connection.close)
        self.assertEqual(connection.execute("PRAGMA foreign_key_check").fetchall(), [])
        self.assertEqual(connection.execute("PRAGMA integrity_check").fetchall(), [("ok",)])

    def test_failed_import_removes_incomplete_output(self):
        (self.data_dir / "parasect_dataset.txt").write_text(
            "domain_id\tsequence\tspecificity\nTest.A1\tAAAA\tmissing substrate\n",
        )
        with self.assertRaisesRegex(ValueError, "No substrate found"):
            create_database(self.output, self.data_dir)
        self.assertFalse(self.output.exists())

    def test_overlapping_source_domains_remain_separate_without_prompts(self):
        with (self.data_dir / "parasect_dataset.txt").open("a") as handle:
            handle.write("Other.A1\tAAA\talanine|valine\n")
        with (self.data_dir / "signatures.fasta").open("a") as handle:
            handle.write(">Other.A1\n" + "A" * 10 + "\n")
        with (self.data_dir / "extended_signatures.fasta").open("a") as handle:
            handle.write(">Other.A1\n" + "A" * 34 + "\n")
        with (self.data_dir / "proteins.fasta").open("a") as handle:
            handle.write(">Other\nMAAAK\n")
        with (self.data_dir / "taxonomy.txt").open("a") as handle:
            handle.write("Other\tBacteria\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tTest\tTest species\tUnknown\n")
        with patch("builtins.input", side_effect=AssertionError("Fresh builds must not prompt")):
            create_database(self.output, self.data_dir)
        engine = create_engine(f"sqlite:///{self.output}")
        self.addCleanup(engine.dispose)
        with Session(engine) as session:
            domains = session.scalars(select(AdenylationDomain)).all()
            self.assertEqual({d.get_name(): [s.name for s in d.substrates] for d in domains}, {
                "Test.A1": ["valine", "alanine", "leucine"],
                "Other.A1": ["alanine", "valine"],
            })

    def test_missing_input_does_not_create_output(self):
        (self.data_dir / "taxonomy.txt").unlink()
        with self.assertRaises(FileNotFoundError):
            create_database(self.output, self.data_dir)
        self.assertFalse(self.output.exists())

    def test_refuses_to_overwrite_an_existing_file(self):
        self.output.write_bytes(b"existing database")
        with self.assertRaises(FileExistsError):
            create_database(self.output, self.data_dir)
        self.assertEqual(self.output.read_bytes(), b"existing database")


if __name__ == "__main__":
    unittest.main()
