"""Substrate order must survive import, correction, and a fresh database session."""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from sqlalchemy import create_engine, select
from sqlalchemy.orm import Session

from parasect.database.build_database import Base, AdenylationDomain, Substrate
from parasect.database.populate_database import create_domain_entries
from parasect.database.process_substrate_corrections import correct_substrate


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


if __name__ == "__main__":
    unittest.main()
