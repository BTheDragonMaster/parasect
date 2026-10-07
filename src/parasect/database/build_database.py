from typing import Any, Optional
from sys import argv

from sqlalchemy import create_engine, CheckConstraint, ForeignKey, JSON
from sqlalchemy.ext.associationproxy import AssociationProxy, association_proxy
from sqlalchemy.ext.orderinglist import ordering_list
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column, relationship


class Base(DeclarativeBase):
    pass


class SubstrateDomainAssociation(Base):
    __tablename__ = "substrate_domain_association"
    __table_args__ = (CheckConstraint("position >= 0"),)

    substrate_name: Mapped[str] = mapped_column(ForeignKey("substrate.name"), primary_key=True)
    domain_id: Mapped[int] = mapped_column(ForeignKey("adenylation_domain.id"), primary_key=True)
    position: Mapped[int]

    substrate: Mapped["Substrate"] = relationship(back_populates="domain_associations")
    domain: Mapped["AdenylationDomain"] = relationship(back_populates="substrate_associations")


substrate_domain_association = SubstrateDomainAssociation.__table__


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

    def __repr__(self):
        return f"{self.protein.get_name()} <-> {self.domain.get_name()}"


class Substrate(Base):
    __tablename__ = "substrate"

    name: Mapped[str] = mapped_column(primary_key=True)
    fingerprint: Mapped[list[int]] = mapped_column(JSON)
    smiles: Mapped[str]
    domain_associations: Mapped[list["SubstrateDomainAssociation"]] = relationship(
        back_populates="substrate", cascade="all, delete-orphan",
    )
    domains: AssociationProxy[list["AdenylationDomain"]] = association_proxy(
        "domain_associations", "domain",
        creator=lambda domain: SubstrateDomainAssociation(domain=domain),
    )

    def __repr__(self):
        return self.name

    def __eq__(self, other):
        if not isinstance(other, Substrate):
            return NotImplemented
        return self.name == other.name

    def __hash__(self):
        return hash(self.name)

    def __lt__(self, other):
        if not isinstance(other, Substrate):
            return NotImplemented
        return self.name < other.name

    def to_json(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "smiles": self.smiles,
            "fingerprint": self.fingerprint
        }


class Taxonomy(Base):
    __tablename__ = "taxonomy"

    id: Mapped[int] = mapped_column(primary_key=True)
    domain: Mapped[str]
    kingdom: Mapped[str]
    phylum: Mapped[str]
    cls: Mapped[str]
    order: Mapped[str]
    family: Mapped[str]
    genus: Mapped[str]
    species: Mapped[str]
    strain: Mapped[Optional[str]] = mapped_column(nullable=True)

    proteins: Mapped[list["Protein"]] = relationship(back_populates="taxonomy")

    def __repr__(self):
        return f"{self.species} {self.strain}"


class Protein(Base):
    __tablename__ = "protein"

    id: Mapped[int] = mapped_column(primary_key=True)

    taxonomy_id: Mapped[int] = mapped_column(ForeignKey("taxonomy.id"))
    taxonomy: Mapped["Taxonomy"] = relationship(back_populates="proteins")

    sequence: Mapped[str]
    synonyms: Mapped[list["ProteinSynonym"]] = relationship(back_populates="protein")
    domains: Mapped[list["ProteinDomainAssociation"]] = relationship(
        back_populates="protein",
        cascade="all, delete-orphan")

    def __repr__(self):
        return self.get_name()

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

    def __repr__(self):
        return self.synonym


class AdenylationDomain(Base):
    __tablename__ = "adenylation_domain"

    id: Mapped[int] = mapped_column(primary_key=True)
    synonyms: Mapped[list["DomainSynonym"]] = relationship(
        back_populates="domain",
        cascade="all, delete-orphan",
        lazy="selectin",
    )
    substrate_associations: Mapped[list["SubstrateDomainAssociation"]] = relationship(
        back_populates="domain",
        cascade="all, delete-orphan",
        order_by="SubstrateDomainAssociation.position",
        collection_class=ordering_list("position"),
    )
    substrates: AssociationProxy[list["Substrate"]] = association_proxy(
        "substrate_associations", "substrate",
        creator=lambda substrate: SubstrateDomainAssociation(substrate=substrate),
    )

    proteins: Mapped[list["ProteinDomainAssociation"]] = relationship(
        back_populates="domain",
        cascade="all, delete-orphan"
    )

    sequence: Mapped[str]
    signature: Mapped[str]
    extended_signature: Mapped[str]

    def __repr__(self):
        return self.get_name()

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

    def __repr__(self):
        return self.synonym

    def to_json(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "synonym": self.synonym,
        }


if __name__ == "__main__":
    engine = create_engine(f"sqlite:///{argv[1]}", echo=True)
    Base.metadata.create_all(engine)
