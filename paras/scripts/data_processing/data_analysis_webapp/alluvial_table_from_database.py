from sys import argv

from sqlalchemy import create_engine
from sqlalchemy.orm import Session

from parasect.database.flatfiles_from_db import taxonomy_to_alluvial


if __name__ == "__main__":
    engine = create_engine(f"sqlite:///{argv[1]}")
    with Session(engine) as session:
        taxonomy_to_alluvial(session, argv[2])
