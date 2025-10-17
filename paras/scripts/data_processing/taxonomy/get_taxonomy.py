from typing import Iterable, TypeVar
from argparse import ArgumentParser, Namespace
import os

from Bio import Entrez
from Bio import SeqIO
from sqlalchemy.orm import Session
from sqlalchemy import create_engine, select

from parasect.database.protein_identifiers import check_id_type, IdType
from parasect.database.build_database import AdenylationDomain, Protein, ProteinSynonym
from parasect.core.writers import write_list


T = TypeVar("T")


def parse_arguments() -> Namespace:
    """Parse arguments from command line

    :return: Arguments
    :rtype: Namespace
    """
    parser = ArgumentParser(description="Retrieve taxonomy for database domains")
    parser.add_argument("-db", '--database', required=True, type=str,
                        help="Path to PARASECT database")
    parser.add_argument("-o", '--out', required=True, type=str,
                        help="Path to output directory")
    parser.add_argument("-e", '--email', required=True, type=str, help="E-mail address")

    arguments = parser.parse_args()

    return arguments


def batch(iterable: Iterable[T], batch_size: int) -> list[list[T]]:
    """Split iterable into batches of size batch_size."""
    lst = list(iterable)
    return [lst[i:i + batch_size] for i in range(0, len(lst), batch_size)]


def get_queryable_proteins(domain: AdenylationDomain) -> list[tuple[Protein, ProteinSynonym]]:
    """Get all associated protein–synonym pairs that are queryable IDs.

    :param domain: adenylation domain
    :type domain: AdenylationDomain
    :return: list of (Protein, ProteinSynonym) pairs
    :rtype: list[tuple[Protein, ProteinSynonym]]
    """
    queryable_pairs: list[tuple[Protein, ProteinSynonym]] = []

    for assoc in domain.proteins:  # ProteinDomainAssociation
        for syn in assoc.protein.synonyms:
            if check_id_type(syn.synonym) in {IdType.REFSEQ, IdType.GENPEPT, IdType.UNIPROT}:
                queryable_pairs.append((assoc.protein, syn))

    return queryable_pairs


def get_best_protein_id(queryable_pairs: list[tuple[Protein, ProteinSynonym]]) -> list[str]:
    """Get best protein synonym for retrieving taxonomy

    :param queryable_pairs: pairs of proteins and protein synonyms
    :type queryable_pairs: list[tuple[Protein, ProteinSynonym]]
    :return: list of protein synonyms that are queryable, with max one synonym representing each protein
    :rtype: str
    """
    priority = {IdType.UNIPROT: 0, IdType.REFSEQ: 1, IdType.GENPEPT: 2}

    # map: protein -> best (protein, synonym) pair
    best_per_protein: dict[Protein, str] = {}

    for protein, synonym in queryable_pairs:
        id_type = check_id_type(synonym.synonym)
        if protein not in best_per_protein:
            best_per_protein[protein] = synonym.synonym
        else:
            # replace only if this synonym has higher priority
            current_best = best_per_protein[protein][1]
            current_priority = priority.get(check_id_type(current_best), 999)
            new_priority = priority.get(id_type, 999)
            if new_priority < current_priority:
                best_per_protein[protein] = synonym.synonym

    # final list of unique proteins with their highest-priority queryable synonym
    unique_queryables = list(best_per_protein.values())
    return unique_queryables


def retrieve_taxonomy(session: Session, out_dir: str, e_mail: str) -> None:
    """
    :param session: PARASECT database session
    :type session: Session
    :param out_dir: path to output directory
    :type out_dir: str
    :param e_mail: e-mail address
    :type e_mail: str
    """
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    out_tax = os.path.join(out_dir, "taxonomy.txt")
    unqueryable_out = os.path.join(out_dir, "unqueryable.txt")
    query = select(AdenylationDomain)
    domains = list(session.scalars(query))
    queryables: list[tuple[Protein, ProteinSynonym]] = []
    unqueryable: list[AdenylationDomain] = []

    for domain in domains:
        queryable_pairs = get_queryable_proteins(domain)
        if not queryable_pairs:
            unqueryable.append(domain)
        else:
            queryables.extend(queryable_pairs)

    write_list([d.get_name() for d in unqueryable], unqueryable_out)

    unique_queryables = get_best_protein_id(queryables)
    unique_queryables.sort()

    batches = batch(unique_queryables, 200)

    Entrez.email = e_mail

    with open(out_tax, 'w') as out:

        for synonym_batch in batches:

            handle = Entrez.efetch(db="protein", id=synonym_batch, rettype="gb", retmode="text")

            for seq_record in SeqIO.parse(handle, "gb"):
                line = f"{seq_record.id}\t" + "\t".join(seq_record.annotations["taxonomy"])
                out.write(f"{line}\n")
                print(line)
        handle.close()


def main():
    args = parse_arguments()

    engine = create_engine(f"sqlite:///{args.database}")

    with Session(engine) as session:
        retrieve_taxonomy(session, args.out, args.email)


if __name__ == "__main__":
    main()
