import os.path
from argparse import ArgumentParser, Namespace
from pprint import pprint

from sqlalchemy.orm import Session
from sqlalchemy import create_engine, select

from parasect.database.build_database import AdenylationDomain, Protein, ProteinDomainAssociation, ProteinSynonym
from parasect.core.featurisation import get_domains
from parasect.core.parsing import parse_fasta_file, write_fasta_file


def parse_arguments() -> Namespace:
    """
    Parse arguments for extracting A domain signatures

    :return: arguments
    :rtype: Namespace
    """
    parser = ArgumentParser(description="Extract signatures for A-domain dataset")
    parser.add_argument('-p', required=True, type=str, help="Path to input protein .fasta")
    parser.add_argument('-t', required=True, type=str, help="Path to temp dir")
    parser.add_argument('-o', required=True, type=str, help="Path to output directory")
    parser.add_argument('-db', required=True, type=str, help="Path to database")

    args = parser.parse_args()

    return args


def extract_signatures(session: Session, protein_file: str, tmp_path: str, out_path: str) -> None:
    """Extract signatures for A-domain dataset

    :param session: database session
    :type session: Session
    :param tmp_path: path to tmp parasect path for storing intermediate files
    :type tmp_path: str
    :param protein_file: Path to full protein sequences
    :type protein_file: str
    :param out_path: path to output directory
    :type out_path: str

    """
    id_to_seq = parse_fasta_file(protein_file)
    original_sequences = set(id_to_seq.keys())

    print("Extracting domains from protein.fasta...")
    domains = get_domains(protein_file, tmp_path, "hmm", "fasta")

    domain_to_db_domain = {}
    db_domain_to_domain = {}
    db_domain_to_positions = {}

    found_db_domains = set()
    all_db_domains = set(session.scalars(select(AdenylationDomain)).all())

    for domain in domains:
        # print(domain.protein_name, domain.domain_nr)
        protein_synonyms = domain.protein_name.split('|')
        # protein_name = protein_synonyms[0]
        #
        # query = (
        #     select(AdenylationDomain, ProteinDomainAssociation)
        #         .join(ProteinDomainAssociation, AdenylationDomain.id == ProteinDomainAssociation.domain_id)
        #         .join(Protein, Protein.id == ProteinDomainAssociation.protein_id)
        #         .join(ProteinSynonym, Protein.id == ProteinSynonym.protein_id)
        #         .where(ProteinSynonym.synonym == protein_name))
        #
        # db_domains = session.execute(query).all()
        match_found = False
        # for db_domain, association in db_domains:
        #     positions = [p for p in domain.extended_signature_positions if p is not None]
        #     sig_start = min(positions) if positions else None
        #     sig_end = max(positions) if positions else None
        #
        #     if sig_start is not None and sig_end is not None:
        #         if association.start - 1 <= sig_start and association.end >= sig_end:
        #             domain_to_db_domain[domain] = db_domain
        #             # if db_domain in db_domain_to_domain:
        #             #     print(f"Multiple domain matches found for {[s.synonym for s in db_domain.synonyms][0]}")
        #             db_domain_to_domain[db_domain] = domain
        #             found_db_domains.add(db_domain)
        #             match_found = True
        #             if domain.extended_signature != db_domain.extended_signature:
        #                 print(f"Mismatching extended signatures for {domain.protein_name} domain {domain.domain_nr}: {domain.extended_signature}, {db_domain.extended_signature}")
        #                 print("Positions:")
        #                 print('\t'.join(map(str, domain.extended_signature_positions)))
        #             if domain.signature != db_domain.signature:
        #                 print(
        #                     f"Mismatching signatures for {domain.protein_name} domain {domain.domain_nr}: {domain.signature}, {db_domain.signature}")
        #
        #             break


        unique_protein_query = select(Protein).join(Protein.synonyms).where(ProteinSynonym.synonym.in_(protein_synonyms)).distinct()
        unique_proteins = session.execute(unique_protein_query).scalars().all()
        for protein in unique_proteins:
            print("protein name", protein.get_name())
            query = (
                select(AdenylationDomain, ProteinDomainAssociation)
                    .join(AdenylationDomain.proteins)
                    .join(ProteinDomainAssociation.protein)
                    .where(Protein.id == protein.id).distinct()
            )

            db_domains = session.execute(query).all()

            for db_domain, association in db_domains:

                name = '|'.join([s.synonym for s in association.protein.synonyms])

                positions = [p for p in domain.extended_signature_positions if p is not None]
                sig_start = min(positions) if positions else None
                sig_end = max(positions) if positions else None

                if sig_start is not None and sig_end is not None:

                    if association.start - 1 <= sig_start and association.end >= sig_end:
                        domain_to_db_domain[domain] = db_domain
                        match_found = True
                        print(domain.protein_name, domain.domain_nr, db_domain.get_name(),
                              association.protein.get_name(), association.id)
                        print("db domain", db_domain.get_name())
                        print("associated protein", association.protein.get_name())
                        print(sig_start, sig_end)
                        if db_domain not in db_domain_to_positions:
                            db_domain_to_positions[db_domain] = {name: []}
                        if name not in db_domain_to_positions[db_domain]:
                            db_domain_to_positions[db_domain][name] = []
                        db_domain_to_positions[db_domain][name].append(domain.extended_signature_positions)

            if not match_found:
                print(f"No database match found for domain {domain.protein_name}.A{domain.domain_nr}")
                print("Positions:")
                print('\t'.join(map(str, domain.extended_signature_positions)))

    missing_domains = all_db_domains - found_db_domains
    # for domain in missing_domains:
    #     print([s.synonym for s in domain.synonyms])
    print(len(all_db_domains))
    print(len(domains))

    pprint(db_domain_to_positions)
    if not os.path.exists(out_path):
        os.mkdir(out_path)

    signatures_out = os.path.join(out_path, 'signatures.fasta')
    extended_out = os.path.join(out_path, 'extended_signatures.fasta')
    positions_file = os.path.join(out_path, 'positions.txt')
    relative_positions_file = os.path.join(out_path, 'relative_positions.txt')
    all_positions_dir = os.path.join(out_path, 'all_positions')
    if not os.path.exists(all_positions_dir):
        os.mkdir(all_positions_dir)

    id_to_sig: dict[str, str] = {}
    id_to_ext: dict[str, str] = {}
    id_to_pos: dict[str, list[int]] = {}

    for db_domain, domain in db_domain_to_domain.items():
        domain_id = '|'.join([s.synonym for s in db_domain.synonyms])
        id_to_sig[domain_id] = domain.signature
        id_to_ext[domain_id] = domain.extended_signature
        id_to_pos[domain_id] = domain.extended_signature_positions

        all_positions_out = os.path.join(all_positions_dir, f"{domain_id}.txt")

        print("Proteins with multiple identical domains")

        for protein_name, positions in db_domain_to_positions[db_domain].items():
            if len(positions) > 1:
                print(protein_name)

    print("Domains with multiple proteins")
    for db_domain, proteins_to_pos in db_domain_to_positions.items():
        if len(proteins_to_pos.keys()) > 1:
            print(db_domain, list[proteins_to_pos.keys()])

        # with open(all_positions_out, 'w') as out:
        #     out.write("protein_name\tdomain_number")
        #     for i in range(34):
        #         out.write(f"\tpos_{i + 1}")
        #     out.write('\n')
        #     for protein_name, positions in db_domain_to_positions[db_domain].items():
        #         out.write(f"{protein_name}")


    write_fasta_file(id_to_sig, signatures_out)
    write_fasta_file(id_to_ext, extended_out)

    with open(positions_file, 'w') as pos_out:
        pos_out.write("domain_id")
        for i in range(34):
            pos_out.write(f"\tpos_{i + 1}")
        pos_out.write("\n")

        domains = sorted(id_to_pos.keys())
        for domain in domains:
            signature_positions = id_to_pos[domain]
            pos_out.write(f"{domain}")
            for i in range(34):
                pos_out.write(f"\t{signature_positions[i]}")

            pos_out.write('\n')


def main() -> None:
    args = parse_arguments()
    engine = create_engine(f"sqlite:///{args.db}")
    try:
        with Session(engine) as session:
            extract_signatures(session, args.p, args.t, args.o)
    finally:
        session.close()


if __name__ == "__main__":
    main()
