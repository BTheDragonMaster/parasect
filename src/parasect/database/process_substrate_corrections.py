from argparse import ArgumentParser, Namespace
import os
import traceback

from sqlalchemy.orm import Session
from sqlalchemy import select, create_engine

from parasect.database.build_database import Substrate
from parasect.database.query_database import get_domains_from_synonym
from parasect.database.populate_database import create_substrate_entries
from parasect.core.parsing import iterate_over_dir, parse_parasect_data


def parse_arguments() -> Namespace:
    """Parse arguments from command line

    :return: User-defined arguments
    :rtype: Namespace
    """
    parser = ArgumentParser(description="Add pending entries to database")
    parser.add_argument('-u', type=str, required=True, help="Path to user submissions folder")
    parser.add_argument('-db', type=str, required=True, help="Path to PARASECT database")

    args = parser.parse_args()
    return args


def correct_substrate(session: Session, domain_synonym: str, correct_substrates: list[str]) -> None:
    """Correct substrate in database

    :param session: database session
    :type session: Session
    :param domain_synonym: one synonym for an adenylation domain
    :type domain_synonym: str
    :param correct_substrates: list of substrates recognised by adenylation domain
    :type correct_substrates: list
    """

    domains = get_domains_from_synonym(session, domain_synonym)
    if not domains:
        raise ValueError(f"No domain found for synonym '{domain_synonym}'")

    domain = domains[0]

    query = select(Substrate).where(Substrate.name.in_(correct_substrates))
    substrates = list(session.scalars(query))

    missing = set(correct_substrates) - {s.name for s in substrates}
    if missing:
        raise ValueError(f"Substrates not found in DB: {', '.join(missing)}")

    domain.substrates.clear()
    domain.substrates.extend(substrates)

    session.add(domain)


def main():
    args = parse_arguments()
    engine = create_engine(f"sqlite:///{args.db}")

    with Session(engine) as session:
        try:
            for _, user_submission_path in iterate_over_dir(args.u, get_dirs=True):
                new_substrates = os.path.join(user_submission_path, "smiles.tsv")
                if os.path.exists(new_substrates):

                    substrates = create_substrate_entries(session, new_substrates)
                    for substrate in substrates:
                        session.add(substrate)

                parasect_path = os.path.join(user_submission_path, "parasect_dataset.txt")
                _, id_to_spec = parse_parasect_data(parasect_path)
                for domain_synonym, substrates in id_to_spec.items():
                    correct_substrate(session, domain_synonym, substrates)

            session.commit()
        except Exception as e:
            print(f"[ERROR] {type(e).__name__}: {e}")
            print("Rolling back changes.")
            traceback.print_exc()
            session.rollback()
        finally:
            session.close()

if __name__ == "__main__":
    main()
