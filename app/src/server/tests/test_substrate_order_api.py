"""The web APIs must retain stored order for a domain's associated substrates."""

import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from flask import Flask
from sqlalchemy import create_engine
from sqlalchemy.orm import Session

from parasect.database.build_database import (
    AdenylationDomain, Base, DomainSynonym, Protein, ProteinDomainAssociation,
    ProteinSynonym, Substrate, Taxonomy,
)

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from routes import compare, network, sql


class TestSubstrateOrderAPI(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.database = Path(temporary.name) / "reference.db"
        self.engine = create_engine(f"sqlite:///{self.database}")
        self.addCleanup(self.engine.dispose)
        Base.metadata.create_all(self.engine)
        self.names = ["valine", "alanine", "leucine"]
        with Session(self.engine) as session:
            substrates = {name: Substrate(name=name, smiles="C", fingerprint=[]) for name in sorted(self.names)}
            taxonomy = Taxonomy(**{
                rank: "Test species" if rank == "species" else "unknown"
                for rank in ["domain", "kingdom", "phylum", "cls", "order", "family", "genus", "species"]
            })
            protein = Protein(sequence="A" * 100, taxonomy=taxonomy,
                              synonyms=[ProteinSynonym(synonym="Test")])
            domain = AdenylationDomain(
                sequence="A" * 50, signature="A" * 10, extended_signature="A" * 34,
                synonyms=[DomainSynonym(synonym="Test.A1")],
                substrates=[substrates[name] for name in self.names],
            )
            session.add(ProteinDomainAssociation(
                protein=protein, domain=domain, domain_number=1, start=0, end=50,
            ))
            session.commit()
            self.domain_id = domain.id

        def get_db():
            with Session(self.engine) as session:
                yield session

        for module in [compare, network]:
            patcher = patch.object(module, "get_db", get_db)
            patcher.start()
            self.addCleanup(patcher.stop)
        for attribute, value in [("_cache", None), ("_cluster_cache", {})]:
            patcher = patch.object(network, attribute, value)
            patcher.start()
            self.addCleanup(patcher.stop)
        patcher = patch.object(sql, "DB_PATH", str(self.database))
        patcher.start()
        self.addCleanup(patcher.stop)

        app = Flask(__name__)
        app.testing = True
        app.register_blueprint(compare.blueprint_compare)
        app.register_blueprint(network.blueprint_network)
        app.register_blueprint(sql.blueprint_sql)
        self.client = app.test_client()

    def test_compare_preserves_order_and_substrate_filter(self):
        response = self.client.post("/api/compare/domains", json={"substrates": ["VALINE"]})
        self.assertEqual(response.status_code, 200)
        rows = response.get_json()["rows"]
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["substrates"], self.names)

    def test_network_cluster_preserves_order_after_session_closes(self):
        response = self.client.get("/api/network/cluster/0")
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.get_json()["nodes"][0]["substrates"], self.names)

    def test_all_dataset_presets_preserve_order(self):
        cases = [
            ("substrate", {"substrate_name": self.names}),
            ("proteinId", {"protein_id": "Test"}),
            ("species", {"species": "Test species"}),
            ("signature", {"signature": "A" * 10, "max_distance": 0}),
        ]
        for preset, params in cases:
            with self.subTest(preset=preset):
                response = self.client.post("/api/sql/preset", json={"preset": preset, "params": params})
                self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
                rows = response.get_json()["rows"]
                self.assertEqual([row["substrate_name"] for row in rows], self.names)


if __name__ == "__main__":
    unittest.main()
