"""Rebuild substrate links with explicit positions, preserving all other data and IDs."""

import csv
import sqlite3
from argparse import ArgumentParser
from collections import defaultdict
from contextlib import closing
from pathlib import Path


def _ordered_associations(connection: sqlite3.Connection, dataset: Path) -> list[tuple[str, int, int]]:
    synonyms = defaultdict(set)
    for synonym, domain_id in connection.execute("SELECT synonym, domain_id FROM domain_synonym"):
        synonyms[synonym].add(domain_id)

    existing = defaultdict(set)
    for name, domain_id in connection.execute(
        "SELECT substrate_name, domain_id FROM substrate_domain_association"
    ):
        existing[domain_id].add(name)

    orders = {}
    with dataset.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"domain_id", "specificity"}.issubset(reader.fieldnames or []):
            raise ValueError("Dataset must contain domain_id and specificity columns")
        for row in reader:
            ids = set().union(*(synonyms[name] for name in row["domain_id"].split("|")))
            if len(ids) != 1:
                raise ValueError(f"Cannot resolve exactly one domain for {row['domain_id']!r}")
            domain_id = ids.pop()
            names = row["specificity"].split("|")
            if len(names) != len(set(names)):
                raise ValueError(f"Duplicate substrates for {row['domain_id']!r}")
            if set(names) != existing[domain_id]:
                raise ValueError(f"Substrates differ from the database for {row['domain_id']!r}")
            if domain_id in orders and orders[domain_id] != names:
                raise ValueError(f"Conflicting substrate orders for domain {domain_id}")
            orders[domain_id] = names

    missing = set(existing) - set(orders)
    if missing:
        raise ValueError(f"Dataset has no substrate order for {len(missing)} database domains")

    return [
        (name, domain_id, position)
        for domain_id, names in orders.items()
        for position, name in enumerate(names)
    ]


def rebuild_substrate_order(database: Path, dataset: Path, output: Path) -> int:
    """Write a new database with source-file substrate order; never modify the input.

    The source must account for every existing association without changing its
    membership. No alphabetical or rowid fallback is used for missing orders.
    Accepts either parasect_dataset.txt or domain_substrate_mapping.txt.
    """
    with closing(sqlite3.connect(database.resolve().as_uri() + "?mode=ro", uri=True)) as source:
        rows = _ordered_associations(source, dataset)
        # Exclusive creation also prevents overwriting the source or an old output.
        with output.open("xb"):
            pass
        try:
            with closing(sqlite3.connect(output)) as destination, destination:
                source.backup(destination)
                destination.execute("PRAGMA foreign_keys = ON")
                destination.execute("BEGIN")
                destination.execute("DROP TABLE substrate_domain_association")
                destination.execute("""
                    CREATE TABLE substrate_domain_association (
                        substrate_name VARCHAR NOT NULL REFERENCES substrate(name),
                        domain_id INTEGER NOT NULL REFERENCES adenylation_domain(id),
                        position INTEGER NOT NULL CHECK (position >= 0),
                        PRIMARY KEY (substrate_name, domain_id)
                    )
                """)
                destination.executemany(
                    "INSERT INTO substrate_domain_association VALUES (?, ?, ?)", rows,
                )
                if destination.execute("PRAGMA foreign_key_check").fetchall():
                    raise ValueError("Rebuilt database failed foreign key validation")
                if destination.execute("PRAGMA integrity_check").fetchall() != [("ok",)]:
                    raise ValueError("Rebuilt database failed integrity validation")
        except Exception:
            output.unlink()
            raise
    return len(rows)


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--database", required=True, type=Path, help="Existing SQLite database")
    parser.add_argument("--dataset", required=True, type=Path, help="Authoritative ordered dataset")
    parser.add_argument("--output", required=True, type=Path, help="New database path (must not exist)")
    args = parser.parse_args()
    count = rebuild_substrate_order(args.database, args.dataset, args.output)
    print(f"Wrote {count} ordered substrate associations to {args.output}")


if __name__ == "__main__":
    main()
