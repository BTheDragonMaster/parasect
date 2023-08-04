import os
import subprocess

from argparse import ArgumentParser

from paras.scripts.parsers.iterate_over_dir import iterate_over_dir


def parse_arguments():
    parser = ArgumentParser()
    parser.add_argument('-i', type=str, required=True, help='Input directory of fasta files.')
    parser.add_argument('-o', type=str, required=True, help="Output directory.")

    args = parser.parse_args()
    return args


def run_fasttree(in_file, out_file):
    with open(out_file, 'w') as out:
        subprocess.run(["FastTree", in_file], stdout=out)


def run():
    args = parse_arguments()
    if not os.path.exists(args.o):
        os.mkdir(args.o)

    for spec, file_path in iterate_over_dir(args.i, '.fasta'):
        out_path = os.path.join(args.o, f"{spec}.nwk")
        run_fasttree(file_path, out_path)


if __name__ == "__main__":
    run()



