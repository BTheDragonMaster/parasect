from typing import Optional, Any
from sys import argv

from sqlalchemy import create_engine, Column, ForeignKey, Table
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column, relationship


class Base(DeclarativeBase):
    pass


substrate_domain_association = Table(
    "substrate_domain_association",
    Base.metadata,
    Column("substrate_name", ForeignKey("substrate.name"), primary_key=True),
    Column("domain_id", ForeignKey("adenylation_domain.id"), primary_key=True),
)


class Substrate(Base):
    __tablename__ = "substrate"

    name: Mapped[str] = mapped_column(primary_key=True)
    smiles: Mapped[str]
    domains: Mapped[list["AdenylationDomain"]] = relationship(secondary=substrate_domain_association,
                                                              back_populates="substrates")

    def to_json(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "smiles": self.smiles
        }


class AdenylationDomain(Base):
    __tablename__ = "adenylation_domain"

    id: Mapped[int] = mapped_column(primary_key=True)
    synonyms: Mapped[list["DomainSynonym"]] = relationship(
        back_populates="domain",
        cascade="all, delete-orphan",
        lazy="selectin",
    )
    ncbi_id: Mapped[Optional[str]]
    substrates: Mapped[list["Substrate"]] = relationship(secondary=substrate_domain_association,
                                                         back_populates="domains")
    sequence: Mapped[str]
    signature: Mapped[str]
    extended_signature: Mapped[str]

    pending: Mapped[bool]

    def to_json(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "ncbi_id": self.ncbi_id,
            "sequence": self.sequence,
            "signature": self.signature,
            "extended_signature": self.extended_signature,
            "pending": self.pending,
            "synonyms": [syn.to_json() for syn in self.synonyms] if self.synonyms else [],
            "substrates": [sub.to_json() for sub in self.substrates] if self.substrates else [],
        }


class DomainSynonym(Base):
    __tablename__ = "domain_synonym"

    id: Mapped[int] = mapped_column(primary_key=True)
    domain_id: Mapped[int] = mapped_column(ForeignKey("adenylation_domain.id"))
    synonym: Mapped[str]

    domain: Mapped[AdenylationDomain] = relationship(back_populates="synonyms")

    def to_json(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "synonym": self.synonym,
        }


if __name__ == "__main__":
    engine = create_engine(f"sqlite:///{argv[1]}", echo=True)
    Base.metadata.create_all(engine)
