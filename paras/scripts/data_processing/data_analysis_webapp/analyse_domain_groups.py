import os
from argparse import ArgumentParser, Namespace

from sqlalchemy.orm import Session

from parasect.database.query_database import get_domains_from_synonym
from parasect.core.parsing import parse_list, iterate_over_dir
from sqlalchemy import create_engine, select
from collections import Counter


def parse_arguments() -> Namespace:
    """Parse arguments from command line

    :return: Arguments
    :rtype: Namespace
    """
    parser = ArgumentParser(description="Retrieve taxonomy for database domains")
    parser.add_argument("-db", '--database', required=True, type=str,
                        help="Path to PARASECT database")
    parser.add_argument("-i", '--input', required=True, type=str,
                        help="Path to directory containing files with domain groups, for instance divided by taxonomic rank")
    parser.add_argument("-o", '--out', required=True, type=str,
                        help="Path to output directory")
    parser.add_argument("-f", "--first_substrate_only", action="store_true",
                        help="If given, store only the first substrate for each domain")

    arguments = parser.parse_args()

    return arguments


def count_sigs_and_specs(session: Session, domain_synonyms: list[str], first_substrate_only: bool) -> tuple[dict[str, int], dict[str, int], dict[str, int]]:
    """Count signatures and specificities within a set of domains

    :param session: PARASECT database session
    :type session: Session
    :param domain_synonyms: list of domain synonyms
    :type domain_synonyms: list[str]
    :param first_substrate_only: if True, count only the first substrate stored for each domain
    :type first_substrate_only: bool
    :return: dictionaries counting specificities and signatures within the domain set
    :rtype: tuple[dict[str, int], dict[str, int], dict[str, int]]
    """

    spec_to_count: dict[str, int] = {}
    sig_to_count: dict[str, int] = {}
    extended_sig_to_count: dict[str, int] = {}
    for synonyms in domain_synonyms:
        synonym = synonyms.split('|')[0]
        domain = get_domains_from_synonym(session, synonym)[0]

        for substrate in domain.substrates:
            if substrate.name not in spec_to_count:
                spec_to_count[substrate.name] = 0
            spec_to_count[substrate.name] += 1

            if first_substrate_only:
                break

        if domain.signature not in sig_to_count:
            sig_to_count[domain.signature] = 0
        if domain.extended_signature not in extended_sig_to_count:
            extended_sig_to_count[domain.extended_signature] = 0

        sig_to_count[domain.signature] += 1
        extended_sig_to_count[domain.extended_signature] += 1

    return spec_to_count, sig_to_count, extended_sig_to_count


def add_to_parent_dict(parent_dict, child_dict, label):
    for key, value in child_dict.items():
        if key in parent_dict:
            parent_dict[key][label] = value
        else:
            parent_dict[key] = {label: value}


def write_parent_dict(parent_dict, out_file, column_names, id_name):

    with open(out_file, 'w') as out:
        out.write(id_name)
        for column_name in column_names:
            out.write(f"\t{column_name}")
        out.write('\n')

        for key in parent_dict:
            out.write(key)
            for column_name in column_names:
                value = 0
                if column_name in parent_dict[key]:
                    value = parent_dict[key][column_name]

                out.write(f"\t{value}")
            out.write('\n')


def write_summary_dict(summary_dict, out_file):
    with open(out_file, 'w') as out:
        for key, value in summary_dict.most_common():
            out.write(f"{key}\t{value}\n")


def main():
    args = parse_arguments()
    if not os.path.exists(args.out):
        os.mkdir(args.out)

    spec_counts_per_group: dict[dict[str, int]] = {}
    sig_counts_per_group: dict[dict[str, int]] = {}
    extended_counts_per_group: dict[dict[str, int]] = {}
    spec_summary: Counter[str] = Counter()
    sig_summary: Counter[str] = Counter()
    extended_summary: Counter[str] = Counter()
    column_names = []

    engine = create_engine(f"sqlite:///{args.database}")
    rank_folder = os.path.join(args.out, 'substrates_per_rank')
    if not os.path.exists(rank_folder):
        os.mkdir(rank_folder)

    with Session(engine) as session:

        for rank_name, file_path in iterate_over_dir(args.input, extension='.txt'):
            column_names.append(rank_name)
            domain_synonyms = parse_list(file_path)
            spec_to_count, sig_to_count, extended_sig_to_count = count_sigs_and_specs(session, domain_synonyms, args.first_substrate_only)

            add_to_parent_dict(spec_counts_per_group, spec_to_count, rank_name)
            add_to_parent_dict(sig_counts_per_group, sig_to_count, rank_name)
            add_to_parent_dict(extended_counts_per_group, extended_sig_to_count, rank_name)

            out_file = os.path.join(rank_folder, f"{rank_name}.txt")
            write_summary_dict(Counter(spec_to_count), out_file)

            spec_summary.update(spec_to_count)
            sig_summary.update(sig_to_count)
            extended_summary.update(extended_sig_to_count)

    column_names.sort()

    spec_path = os.path.join(args.out, 'substrates_per_rank.txt')
    sig_path = os.path.join(args.out, 'signatures_per_rank.txt')
    extended_path = os.path.join(args.out, 'extended_signatures_per_rank.txt')

    spec_summary_path = os.path.join(args.out, 'substrates.txt')
    sig_summary_path = os.path.join(args.out, 'signatures.txt')
    extended_summary_path = os.path.join(args.out, 'extended_signatures.txt')

    write_parent_dict(spec_counts_per_group, spec_path, column_names, "substrate")
    write_parent_dict(sig_counts_per_group, sig_path, column_names, "signature")
    write_parent_dict(extended_counts_per_group, extended_path, column_names, "extended_signature")

    write_summary_dict(spec_summary, spec_summary_path)
    write_summary_dict(sig_summary, sig_summary_path)
    write_summary_dict(extended_summary, extended_summary_path)


if __name__ == "__main__":
    main()
