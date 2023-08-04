from pymol import cmd
import os
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', type=str, required=True, help="Input directory containing homology models of n-terminal portion of A-domain.")
    parser.add_argument('-c', type=str, required=True, help="Input file containing c-terminal domain of 1amu")
    parser.add_argument('-o', type=str, required=True, help="Output directory.")
    args = parser.parse_args()
    return args


def run():
    args = parse_arguments()
    if not os.path.exists(args.o):
        os.mkdir(args.o)
    for file_name in os.listdir(args.n):
        if file_name.endswith('.pdb'):
            file_path = os.path.join(args.n, file_name)
            domain_name = file_name.split('.pdb')[0].replace('|', '___')
            align_to_1amu(file_path, args.c, domain_name, args.o)


def align_to_1amu(pdb_n, pdb_c, domain_name, out_dir):
    n_term = f"{domain_name}_n"
    c_term = f"{domain_name}_c"
    cmd.load(pdb_n, n_term)
    cmd.load(pdb_c, c_term)

    cmd.fetch("1amu")
    cmd.select("chB", "chain B and 1amu")
    cmd.remove("chB")
    cmd.align(n_term, "1amu")
    cmd.align(c_term, "1amu")
    cmd.alter(c_term, "chain='B'")

    cmd.create(domain_name, f"{n_term} | {c_term}")

    file_name = domain_name.replace('___', '|')

    out_path = os.path.join(out_dir, f"{file_name}.pdb")
    cmd.save(out_path, selection=domain_name, format="pdb")
    cmd.reinitialize()


if __name__ == "__main__":
    run()