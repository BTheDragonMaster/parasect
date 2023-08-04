import os
import argparse
from shutil import move
from paras.scripts.parsers.parsers import parse_domain_list


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=str, required=True, help="Path to input pdb folder.")
    parser.add_argument("-d", type=str, required=True, help="Path to domain list.")
    args = parser.parse_args()
    return args


def run():
    args = parse_arguments()
    rename_pdb_and_pse(args.i, args.d)


def rename_pdb_and_pse(pdb_folder, domain_file):
    domain_list = parse_domain_list(domain_file)
    new_name_to_old_name = {}

    for domain in domain_list:
        new_name = domain.replace('|', '_')
        new_name_to_old_name[domain] = new_name

    old_name_to_new_name = {}
    for new_name, old_name in new_name_to_old_name.items():
        old_name_to_new_name[old_name] = new_name

    print(old_name_to_new_name)

    for pdb_file in os.listdir(pdb_folder):
        pdb_path = os.path.join(pdb_folder, pdb_file)
        if os.path.isfile(pdb_path):
            if pdb_file.endswith('.pdb'):
                domain_id = pdb_file.split('.pdb')[0]
                suffix = '.pdb'

            elif pdb_file.endswith('.pse'):
                domain_id = pdb_file.split('.pse')[0]
                suffix = '.pse'
            else:
                continue
            new_name = old_name_to_new_name[domain_id] + suffix
            new_path = os.path.join(pdb_folder, new_name)

            if pdb_path != new_path:
                if '|' not in new_path:
                    print(pdb_path, new_path)
                # move(pdb_path, new_path)


if __name__ == "__main__":
    run()
