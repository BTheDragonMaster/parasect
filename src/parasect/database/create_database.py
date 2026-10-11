"""Create a fresh PARASECT database from the source dataset files."""

import logging
from argparse import ArgumentParser
from pathlib import Path

from sqlalchemy import create_engine
from sqlalchemy.engine import URL
from sqlalchemy.orm import Session

from parasect.database.build_database import Base
from parasect.database.populate_database import populate_db


DEFAULT_DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "database_files"


def create_database(output: Path, data_dir: Path = DEFAULT_DATA_DIR) -> None:
    """Create and populate a new database, preserving each input substrate list's order."""
    files = {
        "parasect_data_path": "parasect_dataset.txt",
        "smiles_path": "smiles.tsv",
        "signature_path": "signatures.fasta",
        "extended_path": "extended_signatures.fasta",
        "protein_path": "proteins.fasta",
        "taxonomy_path": "taxonomy.txt",
    }
    paths = {argument: data_dir / filename for argument, filename in files.items()}
    for path in paths.values():
        if not path.is_file():
            raise FileNotFoundError(f"Missing database input: {path}")

    # A creation command must not silently reuse or overwrite an existing database.
    with output.open("xb"):
        pass
    engine = create_engine(URL.create("sqlite", database=str(output)))
    complete = False
    try:
        with engine.begin() as connection:
            connection.exec_driver_sql("PRAGMA foreign_keys = ON")
            Base.metadata.create_all(connection)
            with Session(bind=connection, autoflush=False) as session:
                # The source dataset already defines distinct domains; do not
                # prompt to merge records merely because their sequences overlap.
                populate_db(session, **{name: str(path) for name, path in paths.items()},
                            check_sequence_overlap=False)
                session.flush()
                if connection.exec_driver_sql("PRAGMA foreign_key_check").fetchall():
                    raise ValueError("Created database failed foreign key validation")
                if connection.exec_driver_sql("PRAGMA integrity_check").fetchall() != [("ok",)]:
                    raise ValueError("Created database failed integrity validation")
        complete = True
    finally:
        engine.dispose()
        if not complete:
            output.unlink()


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path, help="New SQLite database path (must not exist)")
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR,
                        help="Source dataset directory (default: packaged database_files directory)")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    create_database(args.output, args.data_dir)
    print(f"Created database: {args.output}")


if __name__ == "__main__":
    main()
