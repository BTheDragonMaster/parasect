from sqlalchemy import create_engine
from sqlalchemy.orm import Session

from sys import argv
from parasect.core.parsing import parse_list
from parasect.database.query_database import get_domains_from_synonym
from parasect.core.writers import write_list

if __name__ == "__main__":
    engine = create_engine(f"sqlite:///{argv[1]}")
    domain_names = parse_list(argv[2])
    included_substrate_names = parse_list(argv[3])

    promiscuous = []

    with Session(engine) as session:
        domains = []

        for domain_name in domain_names:
            substrates = []
            synonym = domain_name.split('|')[0]
            domain = get_domains_from_synonym(session, synonym)[0]
            for substrate in domain.substrates:
                if substrate.name in included_substrate_names:
                    substrates.append(substrate)
            if len(substrates) > 1:
                promiscuous.append(domain_name)

    write_list(promiscuous, argv[4])


