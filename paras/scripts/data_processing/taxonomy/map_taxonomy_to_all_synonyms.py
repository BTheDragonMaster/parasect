from argparse import ArgumentParser, Namespace

from sqlalchemy.orm import Session
from sqlalchemy import create_engine, select

from parasect.core.parsing import parse_taxonomy_file
from parasect.database.build_database import Protein, ProteinSynonym
from parasect.database.query_database import get_proteins_from_synonym


def parse_arguments() -> Namespace:
    """Parse arguments from command line

    :return: Arguments
    :rtype: Namespace
    """
    parser = ArgumentParser(description="Retrieve taxonomy for database domains")
    parser.add_argument("-db", '--database', required=True, type=str,
                        help="Path to PARASECT database")
    parser.add_argument("-t", '--taxonomy', required=True, type=str,
                        help="Path to taxonomy file")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Path to output file")

    arguments = parser.parse_args()

    return arguments


def rename_proteins(session: Session, taxonomy_file: str, out_file: str) -> None:
    """

    :param session: PARASECT database session
    :type session: Session
    :param taxonomy_file: path to taxonomy file
    :type taxonomy_file: str
    :param out_file: path to output file
    :type out_file: str
    """

    syn_to_tax = parse_taxonomy_file(taxonomy_file)
    synonyms = list(syn_to_tax.keys())

    protein_to_tax = {}

    for synonym in synonyms:
        proteins = get_proteins_from_synonym(session, synonym)
        if proteins:
            protein_name = proteins[0].get_name()
            if protein_name not in protein_to_tax:

                protein_to_tax[protein_name] = syn_to_tax[synonym]
            else:
                print(f"Duplicate proteins: {protein_name}")
        else:
            print(f"No protein found for: {synonym}")

    with open(out_file, 'w') as out:
        out.write(f"protein_id\tdomain\tkingdom\tphylum\tclass\torder\tfamily\tgenus\n")
        proteins = sorted(protein_to_tax.keys())
        for protein in proteins:
            tax = protein_to_tax[protein]
            out.write(f"{protein}\t{tax.domain}\t{tax.kingdom}\t{tax.phylum}\t{tax.cls}\t{tax.order}\t{tax.family}\t{tax.genus}\n")

    subquery = (
        select(ProteinSynonym.protein_id)
            .where(ProteinSynonym.synonym.in_(synonyms))
    )

    # Query proteins for which no synonym has a taxonomy
    proteins = session.scalars(
        select(Protein).where(~Protein.id.in_(subquery))
    ).all()

    for protein in proteins:
        print(protein.get_name())


def main():
    args = parse_arguments()

    engine = create_engine(f"sqlite:///{args.database}")

    with Session(engine) as session:
        rename_proteins(session, args.taxonomy, args.output)


if __name__ == "__main__":
    main()
