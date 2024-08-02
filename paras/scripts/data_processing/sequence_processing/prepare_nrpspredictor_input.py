from sys import argv

from paras.scripts.parsers.fasta import read_fasta


def fasta_to_sigs(fasta_in, sigs_out):
    id_to_seq = read_fasta(fasta_in)
    with open(sigs_out, 'w') as out:
        for seq_id, seq in id_to_seq.items():
            out.write(f"{seq}\t{seq_id}\n")


if __name__ == "__main__":
    fasta_to_sigs(argv[1], argv[2])
