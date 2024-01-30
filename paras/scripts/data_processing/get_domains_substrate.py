import argparse
import os

from paras.scripts.parsers.fasta import read_fasta
from paras.scripts.parsers.parsers import parse_specificities

import paras.data

SPEC_FILE = os.path.join(os.path.dirname(paras.data.__file__), 'trustworthy_parasect_dataset.txt')
SPECIFICITIES = parse_specificities(SPEC_FILE)


def parse_arguments():
    parser = argparse.ArgumentParser(description="Create fasta file containing all domains that recognise a certain substrate.")
    parser.add_argument('-f', type=str, required=True, help="Path to input fasta")
    parser.add_argument('-o', type=str, required=True, help="Path to output fasta")
    parser.add_argument('-s', type=str, required=True, help="Substrate name")

    args = parser.parse_args()
    return args


def run():
    args = parse_arguments()
    id_to_seq = read_fasta(args.f)
    with open(args.o, 'w') as out:
        for seq_id, seq in id_to_seq.items():
            if seq_id in SPECIFICITIES:
                specificities = SPECIFICITIES[seq_id]
                if args.s in specificities:
                    out.write(f">{seq_id}\n{seq}\n")


if __name__ == "__main__":
    run()
