"""Substrate order must survive import, correction, and a fresh database session."""

import sqlite3
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from sqlalchemy import create_engine, select
from sqlalchemy.orm import Session

from parasect.database.build_database import Base, AdenylationDomain, Substrate
from parasect.database.populate_database import create_domain_entries
from parasect.database.process_substrate_corrections import correct_substrate
from parasect.database.rebuild_substrate_order import rebuild_substrate_order


class TestSubstrateOrder(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.engine = create_engine(f"sqlite:///{self.root / 'test.db'}")
        self.addCleanup(self.engine.dispose)
        Base.metadata.create_all(self.engine)
        with Session(self.engine) as session:
            session.add_all([
                Substrate(name=name, smiles="C", fingerprint=[])
                for name in ["alanine", "leucine", "valine"]
            ])
            session.commit()

    def import_domain(self, names):
        dataset = self.root / "dataset.txt"
        signature = self.root / "signature.fasta"
        extended = self.root / "extended.fasta"
        dataset.write_text("domain_id\tsequence\tspecificity\nTest.A1\tAAAA\t" + "|".join(names) + "\n")
        signature.write_text(">Test.A1\nAAAAAAAAAA\n")
        extended.write_text(">Test.A1\n" + "A" * 34 + "\n")
        with Session(self.engine) as session:
            substrates = list(session.scalars(select(Substrate)))
            domains, synonyms = create_domain_entries(
                session, str(dataset), str(signature), str(extended), substrates,
            )
            session.add_all(domains + synonyms)
            session.commit()

    def assert_order(self, expected):
        with Session(self.engine) as session:
            domain = session.scalars(select(AdenylationDomain)).one()
            self.assertEqual([s.name for s in domain.substrates], expected)
            self.assertEqual([a.position for a in domain.substrate_associations], list(range(len(expected))))
            self.assertEqual([s["name"] for s in domain.to_json()["substrates"]], expected)

    def test_import_and_reimport_preserve_source_order(self):
        self.import_domain(["valine", "alanine", "leucine"])
        self.assert_order(["valine", "alanine", "leucine"])
        self.import_domain(["leucine", "valine", "alanine"])
        self.assert_order(["leucine", "valine", "alanine"])

    def test_correction_preserves_requested_order_and_removes_old_links(self):
        self.import_domain(["valine", "alanine"])
        with Session(self.engine) as session:
            correct_substrate(session, "Test.A1", ["leucine", "valine", "leucine"])
            session.commit()
        self.assert_order(["leucine", "valine"])
        with Session(self.engine) as session:
            self.assertEqual(list(session.get(Substrate, "alanine").domains), [])

    def test_insert_and_remove_renumber_positions(self):
        self.import_domain(["valine", "alanine"])
        with Session(self.engine) as session:
            domain = session.scalars(select(AdenylationDomain)).one()
            domain.substrates.insert(1, session.get(Substrate, "leucine"))
            domain.substrates.pop(0)
            session.commit()
        self.assert_order(["leucine", "alanine"])

    def test_shared_substrates_have_independent_domain_orders(self):
        self.import_domain(["valine", "alanine"])
        with Session(self.engine) as session:
            domain = AdenylationDomain(
                sequence="CCCC", signature="C" * 10, extended_signature="C" * 34,
                substrates=[session.get(Substrate, name) for name in ["alanine", "valine"]],
            )
            session.add(domain)
            session.commit()
        with Session(self.engine) as session:
            domains = session.scalars(select(AdenylationDomain).order_by(AdenylationDomain.id)).all()
            self.assertEqual([[s.name for s in d.substrates] for d in domains],
                             [["valine", "alanine"], ["alanine", "valine"]])


class TestRebuildSubstrateOrder(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.database = self.root / "original.db"
        self.output = self.root / "ordered.db"
        self.dataset = self.root / "dataset.txt"
        connection = sqlite3.connect(self.database)
        self.addCleanup(connection.close)
        connection.executescript("""
            CREATE TABLE adenylation_domain (id INTEGER PRIMARY KEY, sequence TEXT);
            CREATE TABLE substrate (name TEXT PRIMARY KEY);
            CREATE TABLE domain_synonym (synonym TEXT, domain_id INTEGER);
            CREATE TABLE substrate_domain_association (
                substrate_name TEXT REFERENCES substrate(name),
                domain_id INTEGER REFERENCES adenylation_domain(id),
                PRIMARY KEY (substrate_name, domain_id)
            );
            INSERT INTO adenylation_domain VALUES (42, 'AAAA'), (57, 'CCCC');
            INSERT INTO substrate VALUES ('alanine'), ('leucine'), ('valine');
            INSERT INTO domain_synonym VALUES ('Test.A1', 42), ('Alias.A1', 42), ('Other.A1', 57);
            INSERT INTO substrate_domain_association VALUES
                ('alanine', 42), ('leucine', 42), ('valine', 42), ('alanine', 57);
        """)
        self.source_bytes = self.database.read_bytes()
        self.write_dataset()

    def write_dataset(self, first="valine|alanine|leucine", second="Other.A1\talanine\n"):
        self.dataset.write_text(
            "domain_id\tspecificity\n" + f"Test.A1|Alias.A1\t{first}\n" + second,
        )

    def test_rebuild_preserves_ids_and_data_and_leaves_source_untouched(self):
        self.assertEqual(rebuild_substrate_order(self.database, self.dataset, self.output), 4)
        connection = sqlite3.connect(self.output)
        self.addCleanup(connection.close)
        self.assertEqual(connection.execute(
            "SELECT substrate_name, position FROM substrate_domain_association "
            "WHERE domain_id=42 ORDER BY position"
        ).fetchall(), [("valine", 0), ("alanine", 1), ("leucine", 2)])
        self.assertEqual(connection.execute("SELECT * FROM adenylation_domain ORDER BY id").fetchall(),
                         [(42, "AAAA"), (57, "CCCC")])
        self.assertEqual(connection.execute("PRAGMA foreign_key_check").fetchall(), [])
        self.assertEqual(self.database.read_bytes(), self.source_bytes)

    def test_rebuild_rejects_missing_mismatched_duplicate_and_conflicting_orders(self):
        cases = [
            ("valine|alanine", "Other.A1\talanine\n"),
            ("valine|alanine|leucine", ""),
            ("valine|alanine|leucine|valine", "Other.A1\talanine\n"),
            ("valine|alanine|leucine", "Alias.A1\talanine|leucine|valine\nOther.A1\talanine\n"),
            ("valine|alanine|leucine", "Unknown.A1\talanine\n"),
        ]
        for first, second in cases:
            with self.subTest(first=first, second=second):
                self.write_dataset(first, second)
                with self.assertRaises(ValueError):
                    rebuild_substrate_order(self.database, self.dataset, self.output)
                self.assertFalse(self.output.exists())
                self.assertEqual(self.database.read_bytes(), self.source_bytes)

    def test_rebuild_refuses_to_overwrite_input(self):
        with self.assertRaises(FileExistsError):
            rebuild_substrate_order(self.database, self.dataset, self.database)
        self.assertEqual(self.database.read_bytes(), self.source_bytes)


if __name__ == "__main__":
    unittest.main()
