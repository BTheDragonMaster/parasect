from sys import argv
from paras.scripts.parsers.fasta import read_fasta


def sort_sequences(fasta_file, out_file):
    id_to_seq = read_fasta(fasta_file)
    with open(out_file, 'w') as out:
        for seq_id in sorted(id_to_seq.keys()):
            seq = id_to_seq[seq_id]
            out.write(f">{seq_id}\n{seq}\n")


if __name__ == "__main__":
    sort_sequences(argv[1], argv[2])