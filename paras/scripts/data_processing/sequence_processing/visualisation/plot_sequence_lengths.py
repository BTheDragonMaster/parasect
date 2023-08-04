import matplotlib.pyplot as plt
from sys import argv
from paras.scripts.parsers.fasta import read_fasta


def plot_seqlengths(fasta_file):
    seqlengths = []
    id_to_seq = read_fasta(fasta_file)
    for sequence in id_to_seq.values():
        seqlengths.append(len(sequence))

    plt.hist(seqlengths, 20)
    plt.show()


if __name__ == "__main__":
    plot_seqlengths(argv[1])
    
