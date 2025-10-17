from sys import argv

from sqlalchemy import create_engine
from sqlalchemy.orm import Session

from parasect.database.flatfiles_from_db import substrates_to_cytoscape

if __name__ == "__main__":
    engine = create_engine(f"sqlite:///{argv[1]}")

    with Session(engine) as session:
        substrates_to_cytoscape(session, argv[2])