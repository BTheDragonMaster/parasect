from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from .constants import DB_PATH

# Absolute path: the URL used to be a bare "sqlite:///parasect.db", resolved
# against the process's working directory, which only worked because a second
# copy of the database was baked into /app alongside the packaged one.
DATABASE_URL = f"sqlite:///{DB_PATH}"

engine = create_engine(
    DATABASE_URL,
    connect_args={"check_same_thread": False},  # for SQLite only
    pool_size=10,
    max_overflow=20,
    pool_timeout=30,
)

SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)


def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()
