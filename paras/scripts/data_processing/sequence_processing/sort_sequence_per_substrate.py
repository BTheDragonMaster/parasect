import os
import argparse

from paras.scripts.parsers.parsers import parse_domain_list
from paras.scripts.parsers.fasta import read_fasta, write_fasta


def parse_arguments():
    parser = argparse.ArgumentParser(description="Script to make separate datasets for each substrate.")
    parser.add_argument('-d', type=str, required=True, help='Directory to domain list stored per substrate.')
    parser.add_argument('-f', type=str, required=True, help="Fasta file.")
    parser.add_argument('-o', type=str, required=True, help="Output directory.")

    args = parser.parse_args()
    return args


def remove_gaps(sequences):
    positions = []
    for i, char in enumerate(sequences[0]):
        for sequence in sequences:
            if sequence[i] != '-':
                positions.append(i)
                break

    new_sequences = []

    for sequence in sequences:
        new_sequence = []
        for position in positions:
            new_sequence.append(sequence[position])
        new_sequence = ''.join(new_sequence)
        new_sequences.append(new_sequence)

    print(len(sequences[0]), len(positions))

    return new_sequences


def run():
    args = parse_arguments()
    domain_to_seq = read_fasta(args.f)
    if not os.path.exists(args.o):
        os.mkdir(args.o)
    for file_name in os.listdir(args.d):
        if file_name.endswith('.txt'):
            spec = file_name.split('.txt')[0]
            file_path = os.path.join(args.d, file_name)
            domain_list = parse_domain_list(file_path)
            sequences = []
            for domain in domain_list:
                sequences.append(domain_to_seq[domain])

            sequences = remove_gaps(sequences)
            new_domain_to_seq = dict(zip(domain_list, sequences))
            out_path = os.path.join(args.o, f"{spec}.fasta")
            write_fasta(new_domain_to_seq, out_path)


if __name__ == "__main__":
    run()