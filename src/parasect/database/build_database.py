from typing import Any
from sys import argv

from sqlalchemy import create_engine, Column, ForeignKey, Table, String, JSON
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column, relationship


class Base(DeclarativeBase):
    pass


substrate_domain_association = Table(
    "substrate_domain_association",
    Base.metadata,
    Column("substrate_name", String, ForeignKey("substrate.name"), primary_key=True),
    Column("domain_id", ForeignKey("adenylation_domain.id"), primary_key=True),
)


class ProteinDomainAssociation(Base):
    __tablename__ = "protein_domain_association"

    id: Mapped[int] = mapped_column(primary_key=True)
    protein_id: Mapped[int] = mapped_column(ForeignKey("protein.id"))
    domain_id: Mapped[int] = mapped_column(ForeignKey("adenylation_domain.id"))
    domain_number: Mapped[int]

    start: Mapped[int]
    end: Mapped[int]

    protein: Mapped["Protein"] = relationship(back_populates="domains")
    domain: Mapped["AdenylationDomain"] = relationship(back_populates="proteins")


class Substrate(Base):
    __tablename__ = "substrate"

    name: Mapped[str] = mapped_column(primary_key=True)
    fingerprint: Mapped[list[int]] = mapped_column(JSON)
    smiles: Mapped[str]
    domains: Mapped[list["AdenylationDomain"]] = relationship(secondary=substrate_domain_association,
                                                              back_populates="substrates")

    def to_json(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "smiles": self.smiles,
            "fingerprint": self.fingerprint
        }


class Protein(Base):
    __tablename__ = "protein"

    id: Mapped[int] = mapped_column(primary_key=True)
    sequence: Mapped[str]
    synonyms: Mapped[list["ProteinSynonym"]] = relationship(back_populates="protein")
    domains: Mapped[list["ProteinDomainAssociation"]] = relationship(
        back_populates="protein",
        cascade="all, delete-orphan")

    def get_name(self) -> str:
        """
        Return the protein name
        :return: protein name
        :rtype: str
        """
        return '|'.join([s.synonym for s in self.synonyms])


class ProteinSynonym(Base):
    __tablename__ = "protein_synonym"

    id: Mapped[int] = mapped_column(primary_key=True)
    protein_id: Mapped[int] = mapped_column(ForeignKey("protein.id"))
    synonym: Mapped[str]

    protein: Mapped[Protein] = relationship(back_populates="synonyms")


class AdenylationDomain(Base):
    __tablename__ = "adenylation_domain"

    id: Mapped[int] = mapped_column(primary_key=True)
    synonyms: Mapped[list["DomainSynonym"]] = relationship(
        back_populates="domain",
        cascade="all, delete-orphan",
        lazy="selectin",
    )
    substrates: Mapped[list["Substrate"]] = relationship(secondary=substrate_domain_association,
                                                         back_populates="domains")

    proteins: Mapped[list["ProteinDomainAssociation"]] = relationship(
        back_populates="domain",
        cascade="all, delete-orphan"
    )

    sequence: Mapped[str]
    signature: Mapped[str]
    extended_signature: Mapped[str]

    def get_name(self) -> str:
        """
        Return the domain name
        :return: domain name
        :rtype: str
        """
        return '|'.join([s.synonym for s in self.synonyms])

    def to_json(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "sequence": self.sequence,
            "signature": self.signature,
            "extended_signature": self.extended_signature,
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
