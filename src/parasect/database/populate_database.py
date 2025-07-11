from argparse import ArgumentParser, Namespace
from sqlalchemy import create_engine
from sqlalchemy.orm import Session

from parasect.core.parsing import load_parasect_data
from parasect.database.build_database import AdenylationDomain, Substrate, DomainSynonym


def parse_args() -> Namespace:
    parser = ArgumentParser(description="Populate empty SQL database with data")
    parser.add_argument("--parasect", required=True, type=str, help="Path to parasect dataset file")
    parser.add_argument("--smiles", required=True, type=str, help="Path to parasect SMILES file")
    parser.add_argument("--database", required=True, type=str, help="Path to parasect database")
    parser.add_argument("--signature", required=True, type=str, help="Path to fasta containing A-domain signatures")
    parser.add_argument("--extended", required=True, type=str,
                        help="Path to fasta containing A-domain extended signatures")
    args = parser.parse_args()

    return args


def populate_db(database_path: str, parasect_data_path: str, smiles_path: str, signature_path: str, extended_path: str):
    engine = create_engine(f"sqlite:///{database_path}")
    domains, substrates = load_parasect_data(parasect_data_path, smiles_path, signature_path, extended_path)

    with Session(engine) as session:
        domain_entries = []
        substrate_entries = []
        domain_synonyms = []

        for substrate in substrates:
            substrate_entries.append(Substrate(name=substrate.name, smiles=substrate.smiles))

        for domain in domains:
            domain_entry = AdenylationDomain(id=domain.id, ncbi_id=None,
                                             sequence=domain.sequence, signature=domain.signature,
                                             extended_signature=domain.extended_signature,
                                             pending=False)
            for synonym in domain.synonyms:
                domain_synonyms.append(DomainSynonym(synonym=synonym, domain=domain_entry))

            for substrate_entry in substrate_entries:
                if substrate_entry.name in [s.name for s in domain.substrates]:
                    domain_entry.substrates.append(substrate_entry)

            domain_entries.append(domain_entry)

        entries = domain_entries + substrate_entries + domain_synonyms
        session.add_all(entries)
        session.commit()


def main():
    args = parse_args()
    populate_db(args.database, args.parasect, args.smiles, args.signature, args.extended)


if __name__ == "__main__":
    main()
