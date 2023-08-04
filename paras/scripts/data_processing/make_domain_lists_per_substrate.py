import os
import argparse

from paras.scripts.parsers.parsers import parse_specificities, parse_domain_list
from paras.scripts.data_analysis.count_substrates import count_substrates_from_file
from paras.scripts.data_processing.relabel_data import relabel_from_inclusion_limit



def parse_arguments():
    parser = argparse.ArgumentParser(description="Script to make separate datasets for each substrate.")
    parser.add_argument('-s', type=str, required=True, help='Path to file containing specificity data.')
    parser.add_argument('-d', type=str, required=True, help='Path to file containing domain list.')
    parser.add_argument('-o', type=str, required=True, help="Output directory.")
    parser.add_argument('-l', type=int, required=True, help="Inclusion limit.")

    args = parser.parse_args()
    return args


def run():
    args = parse_arguments()
    domain_to_spec = parse_specificities(args.s)
    domain_list = parse_domain_list(args.d)
    substrate_counts = count_substrates_from_file(args.d, args.s)
    domain_to_filtered = relabel_from_inclusion_limit(domain_list, domain_to_spec, substrate_counts, args.l)

    spec_to_domains = {}
    for domain, specs in domain_to_filtered.items():
        for spec in specs:
            if spec not in spec_to_domains:
                spec_to_domains[spec] = []
            spec_to_domains[spec].append(domain)

    if not os.path.exists(args.o):
        os.mkdir(args.o)

    for spec, domains in spec_to_domains.items():
        out_file = os.path.join(args.o, f"{spec}.txt")
        with open(out_file, 'w') as out:
            for domain in domains:
                out.write(f"{domain}\n")


if __name__ == "__main__":
    run()