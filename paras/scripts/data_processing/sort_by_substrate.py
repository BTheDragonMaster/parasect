import os
import argparse

from paras.scripts.parsers.parsers import parse_specificities


def parse_arguments():
    parser = argparse.ArgumentParser(description="Script to make separate datasets for each substrate.")
    parser.add_argument('-i', type=str, required=True, help='Path to file containing specificity data.')
    parser.add_argument('-o', type=str, required=True, help="Output directory.")

    args = parser.parse_args()
    return args


def run():
    args = parse_arguments()
    sort_by_substrate(args.i, args.o)


def sort_by_substrate(specificities_file, out_dir):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    domain_to_specificities = parse_specificities(specificities_file)
    specificity_to_domain = {}
    all_specificities = set()
    for domain, specificities in domain_to_specificities.items():
        for specificity in specificities:
            if specificity not in specificity_to_domain:
                specificity_to_domain[specificity] = []
            specificity_to_domain[specificity].append(domain)
            all_specificities.add(specificity)

    for specificity in all_specificities:
        file_name = os.path.join(out_dir, f"{specificity}.txt")
        with open(file_name, 'w') as out_file:
            out_file.write("domain\tactivates_substrate\n")
            for domain, specificities in domain_to_specificities.items():
                if specificity in specificities:
                    out_file.write(f"{domain}\t1\n")
                else:
                    out_file.write(f"{domain}\t0\n")


if __name__ == "__main__":
    run()
    