import os
from argparse import ArgumentParser, Namespace
from collections import Counter

from sqlalchemy import create_engine, select
from sqlalchemy.orm import Session


from parasect.core.chem import get_jaccard_distance_from_fingerprint
from parasect.database.build_database import AdenylationDomain, Protein, Substrate
from parasect.database.query_database import get_domains_from_synonym
from parasect.core.writers import write_list
from parasect.core.parsing import parse_list


def parse_arguments() -> Namespace:
    """ Parse command line arguments
    :return: arguments
    :rtype: Namespace
    """
    parser = ArgumentParser(description="Retrieve flatfiles from database")
    parser.add_argument('-db', required=True, type=str, help="Path to database")
    parser.add_argument('-o', required=True, type=str, help="Path to output directory")

    args = parser.parse_args()
    return args


def substrates_to_cytoscape(session: Session, out_dir: str) -> None:
    """Write substrate tables for cytoscape analysis from database

    :param session: database session
    :type session: Session
    :param out_dir: output directory
    :type out_dir: str
    """
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    counts_out = os.path.join(out_dir, "counts.txt")
    distances_out = os.path.join(out_dir, "distances.txt")
    all_substrates = list(session.scalars(select(Substrate)).all())

    with open(counts_out, 'w') as out:

        for substrate in all_substrates:
            count = len(substrate.domains)

            if count > 100:
                out.write(f"{substrate.name}\t4\n")
            elif count > 50:
                out.write(f"{substrate.name}\t3\n")
            elif count > 10:
                out.write(f"{substrate.name}\t2\n")
            elif count > 0:
                out.write(f"{substrate.name}\t1\n")
    with open(distances_out, 'w') as out:
        out.write("substrate_1\tsubstrate_2\tdistance\n")
        for i, substrate_1 in enumerate(all_substrates):
            for j, substrate_2 in enumerate(all_substrates):
                if j > i and len(substrate_1.domains) > 0 and len(substrate_2.domains) > 0:
                    distance = get_jaccard_distance_from_fingerprint(set(substrate_1.fingerprint), set(substrate_2.fingerprint))
                    out.write(f"{substrate_1.name}\t{substrate_2.name}\t{distance}\n")


def taxonomy_to_alluvial(session: Session, out_file: str) -> None:
    """
    Write taxonomy distribution to file which counts occurrences of domain-kingdom-phylum combinations

    :param session: database session
    :type session: Session
    :param out_file: path to output file
    :type out_file: str
    """
    tax_list = []
    for domain in list(session.scalars(select(AdenylationDomain)).all()):
        taxonomy = domain.proteins[0].protein.taxonomy
        domain = taxonomy.domain if taxonomy.domain is not None else "Unknown"
        kingdom = taxonomy.kingdom if taxonomy.kingdom is not None else "Unknown"
        phylum = taxonomy.phylum if taxonomy.phylum is not None else "Unknown"
        tax_list.append((domain, kingdom, phylum))

    counter = Counter(tax_list)
    with open(out_file, 'w') as out:
        out.write(f"Number\tDomain\tKingdom\tPhylum\n")
        for tax, count in counter.items():
            out.write(f"{count}\t{tax[0]}\t{tax[1]}\t{tax[2]}\n")


def write_substrate_names(session: Session, out_path: str) -> None:
    """Write substrates from database to file

    :param session: database session
    :type session: Session
    :param out_path: path to output directory
    :type out_path: str
    """
    all_substrates = list(session.scalars(select(Substrate)).all())
    substrate_names = []
    for substrate in all_substrates:
        if len(substrate.domains) > 0:
            substrate_names.append(substrate.name)
    write_list(substrate_names, out_path)


def write_domain_names(session: Session, out_path: str) -> None:
    all_domains = list(session.scalars(select(AdenylationDomain)).all())
    write_list([d.get_name() for d in all_domains], out_path)


def write_domains_fasta(session: Session, domain_path: str, out_path: str) -> None:
    domain_names = parse_list(domain_path)
    domains = []
    for name in domain_names:
        domains.append(get_domains_from_synonym(session, name.split('|')[0])[0])

    with open(out_path, 'w') as out:
        for domain in domains:
            out.write(f">{domain.get_name()}\n{domain.sequence}\n")


def db_to_flatfiles(database_path: str, out_dir: str) -> None:
    """Write database entries to machine-readable flatfiles (.fasta, .tsv, .txt)

    :param database_path: path to sqlite database
    :type database_path: str
    :param out_dir: path to output directory
    :type out_dir: str

    """

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    engine = create_engine(f"sqlite:///{database_path}")
    domains_path = os.path.join(out_dir, "domains.fasta")
    proteins_path = os.path.join(out_dir, "proteins.fasta")
    substrates_path = os.path.join(out_dir, "smiles.tsv")
    parasect_path = os.path.join(out_dir, "parasect_dataset.txt")
    signatures_path = os.path.join(out_dir, "signatures.fasta")
    extended_signatures_path = os.path.join(out_dir, "extended_signatures.fasta")

    with Session(engine) as session:
        domains: list[AdenylationDomain] = list(session.scalars(select(AdenylationDomain)).all())
        proteins: list[Protein] = list(session.scalars(select(Protein)).all())
        substrates: list[Substrate] = list(session.scalars(select(Substrate)).all())

        domains.sort(key=lambda x: x.synonyms[0].synonym)
        proteins.sort(key=lambda x: x.synonyms[0].synonym)
        substrates.sort(key=lambda x: x.name)

        with open(domains_path, 'w') as domains_out, open(parasect_path, 'w') as parasect_out,\
                open(signatures_path, 'w') as signatures_out, open(extended_signatures_path, 'w') as ext_out:
            parasect_out.write("domain_id\tsequence\tspecificity\n")
            for domain in domains:
                domain_id = '|'.join([s.synonym for s in domain.synonyms])
                substrate_name = '|'.join([s.name for s in domain.substrates])
                domains_out.write(f">{domain_id}\n{domain.sequence}\n")
                parasect_out.write(f"{domain_id}\t{domain.sequence}\t{substrate_name}\n")
                signatures_out.write(f">{domain_id}\n{domain.signature}\n")
                ext_out.write(f">{domain_id}\n{domain.extended_signature}\n")
        with open(proteins_path, 'w') as proteins_out:
            for protein in proteins:
                protein_id = '|'.join([s.synonym for s in protein.synonyms])
                proteins_out.write(f">{protein_id}\n{protein.sequence}\n")
        with open(substrates_path, 'w') as substrates_out:
            substrates_out.write("substrate\tsmiles\n")
            for substrate in substrates:
                if len(substrate.domains) > 0:
                    substrates_out.write(f"{substrate.name}\t{substrate.smiles}\n")


if __name__ == "__main__":
    args = parse_arguments()
    db_to_flatfiles(args.db, args.o)

