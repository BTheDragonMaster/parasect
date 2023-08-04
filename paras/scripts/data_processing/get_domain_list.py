import os
import argparse

from paras.scripts.parsers.tabular import Tabular
from paras.scripts.parsers.parsers import parse_domain_list


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=str, required=True, help="Path to input file or pdb folder.")
    parser.add_argument("-o", type=str, required=True, help="Path to output file.")
    parser.add_argument("-d", type=str, default=None, help="Path to file containing reference IDs")
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-folder', action='store_true', help="Specifies input as a folder")
    group.add_argument('-file', action='store_true', help="Specifies input as a file")
    args = parser.parse_args()
    return args


def run():
    args = parse_arguments()
    if args.folder:
        from_folder(args.i, args.o, args.d)
    elif args.file:
        from_file(args.i, args.o)
    else:
        raise Exception("Specify if input is a file or a folder.")


def from_file(parasect_file, out_file):
    parasect_data = Tabular(parasect_file, [0])
    domain_ids = parasect_data.get_column("domain_id")
    write_domain_ids(domain_ids, out_file)


def from_folder(pdb_folder, out_file, ref_domain_file):
    ref_to_folder = {}
    for domain in parse_domain_list(ref_domain_file):
        new_name = domain.replace('|', '_')
        ref_to_folder[domain] = new_name

    folder_to_ref = {}
    for ref, folder in ref_to_folder.items():
        folder_to_ref[folder] = ref
    domain_ids = []
    for pdb_file in os.listdir(pdb_folder):
        pdb_path = os.path.join(pdb_folder, pdb_file)
        if pdb_file.endswith('.pdb') and os.path.isfile(pdb_path):
            domain_id = pdb_file.split('.pdb')[0]
            domain_ids.append(folder_to_ref[domain_id])
    write_domain_ids(domain_ids, out_file)


def make_domain_mapping(ref_domain_file):
    domain_list = parse_domain_list(ref_domain_file)
    new_name_to_old_name = {}

    for domain in domain_list:
        new_name = domain.replace('|', '_')
        new_name_to_old_name[domain] = new_name

    old_name_to_new_name = {}
    for new_name, old_name in new_name_to_old_name.items():
        old_name_to_new_name[old_name] = new_name

    return old_name_to_new_name


def write_domain_ids(domain_ids, out_file):
    with open(out_file, 'w') as out:
        for domain_id in domain_ids:
            out.write(f"{domain_id}\n")


if __name__ == "__main__":
    run()
