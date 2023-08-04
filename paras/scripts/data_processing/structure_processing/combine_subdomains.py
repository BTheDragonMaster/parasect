from pymol import cmd
import os
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', type=str, required=True, help="Input directory containing homology models of n-terminal portion of A-domain.")
    parser.add_argument('-c', type=str, required=True, help="Input directory containing of c-terminal portion of A-domain.")
    parser.add_argument('-o', type=str, required=True, help="Output directory.")
    args = parser.parse_args()
    return args


def run():
    args = parse_arguments()
    if not os.path.exists(args.o):
        os.mkdir(args.o)
    pdb_path_pairs = make_file_pairs(args.n, args.c)
    for domain_name, paths in pdb_path_pairs.items():
        pdb_n, pdb_c = paths
        align_to_1amu(pdb_n, pdb_c, domain_name, args.o)


def make_file_pairs(n_term_dir, c_term_dir):
    pdb_path_pairs = {}
    for pdb_file in os.listdir(c_term_dir):
        if pdb_file.endswith('.pdb'):
            # domain_name = pdb_file.split('_0.pdb')[0]
            domain_name = pdb_file.split('.pdb')[0]
            domain_name = domain_name.replace('|', '___')

            pdb_path_c = os.path.join(c_term_dir, pdb_file)
            pdb_path_n = os.path.join(n_term_dir, pdb_file)

            print(pdb_path_c)
            print(pdb_path_n)

            assert os.path.isfile(pdb_path_n)
            assert os.path.isfile(pdb_path_c)

            pdb_path_pairs[domain_name] = (pdb_path_n, pdb_path_c)

    return pdb_path_pairs


def align_to_1amu(pdb_n, pdb_c, domain_name, out_dir):
    n_term = f"{domain_name}_n"
    c_term = f"{domain_name}_c"
    cmd.load(pdb_n, n_term)
    cmd.load(pdb_c, c_term)

    cmd.fetch("1amu")
    cmd.select("chB", "chain B")
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