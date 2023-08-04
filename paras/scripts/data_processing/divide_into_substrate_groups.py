import os
import argparse

from paras.scripts.parsers.parsers import parse_specificities, parse_domain_list
from paras.scripts.data_analysis.count_substrates import count_substrates_from_file
from paras.scripts.data_processing.relabel_data import relabel_from_inclusion_limit

SUBSTRATE_GROUPS = {"small_hydrophobic": ["isoleucine", "leucine", "valine"],
                    "aromatic": ["tryptophan", "tyrosine", "phenylalanine", "histidine", "R-beta-hydroxytyrosine"],
                    "aromatic_acids": ["anthranilic acid", "salicylic acid", "2,3-dihydroxybenzoic acid"],
                    "phenylglycine": ["4-hydroxyphenylglycine", "3,5-dihydroxyphenylglycine"],
                    "asx_glx": ["asparagine", "aspartic acid", "glutamine", "glutamic acid", "2-aminoadipic acid"],
                    "cyclic_aliphatic": ["proline", "pipecolic acid"],
                    "NH2_containing": ["ornithine", "lysine", "arginine", "N5-hydroxyornithine",
                                       "N5-formyl-N5-hydroxyornithine"],
                    "small_polar": ["cysteine", "threonine", "serine", "2,4-diaminobutyric acid"],
                    "small": ["glycine", "alanine", "beta-alanine", "D-alanine", "2-aminoisobutyric acid"]}


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

    group_to_domains = {}
    for domain, specs in domain_to_filtered.items():
        for spec in specs:
            for group, substrates in SUBSTRATE_GROUPS.items():
                if spec in substrates:
                    if group not in group_to_domains:
                        group_to_domains[group] = []
                    if domain not in group_to_domains[group]:
                        group_to_domains[group].append(domain)

    if not os.path.exists(args.o):
        os.mkdir(args.o)

    for group, domains in group_to_domains.items():
        out_file = os.path.join(args.o, f"{group}.txt")
        with open(out_file, 'w') as out:
            for domain in domains:
                out.write(f"{domain}\n")


if __name__ == "__main__":
    run()