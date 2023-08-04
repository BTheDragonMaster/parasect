import os
import argparse
import shutil

from paras.scripts.parsers.parsers import parse_paras_specificities


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, required=True, help="Input directory containing pdb files of A-domains.")
    parser.add_argument('-o', type=str, required=True, help="Output directory.")
    parser.add_argument('-s', type=str, required=True, help="Path to specificities file.")
    args = parser.parse_args()
    return args


def run():
    args = parse_arguments()
    if not os.path.exists(args.o):
        os.mkdir(args.o)

    domain_to_specificities = parse_paras_specificities(args.s)
    pdb_name_to_specificities = {}
    pdb_name_to_domain = {}

    for domain, specs in domain_to_specificities.items():
        pdb_name = domain.replace('|', '_')
        pdb_name_to_specificities[pdb_name] = specs
        pdb_name_to_domain[pdb_name] = domain

    for pdb_file in os.listdir(args.i):
        pdb_path = os.path.join(args.i, pdb_file)
        if os.path.isfile(pdb_path) and pdb_file.endswith('.pdb'):
            if '_out.pdb' in pdb_file:
                domain_name = pdb_file.split('_out.pdb')[0]
            else:
                domain_name = pdb_file.split('.pdb')[0]

            specificity = pdb_name_to_specificities[domain_name][0]
            folder_path = os.path.join(args.o, specificity)
            if not os.path.exists(folder_path):
                os.mkdir(folder_path)
            domain = pdb_name_to_domain[domain_name]

            out_path = os.path.join(folder_path, f"{domain}.pdb")
            shutil.copyfile(pdb_path, out_path)


if __name__ == "__main__":
    run()
