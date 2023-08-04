import matplotlib.pyplot as plt
from sys import argv


def plot_bitscores(network_file):
    bitscores = []
    with open(network_file, 'r') as network:
        network.readline()
        for line in network:
            info = line.split('\t')
            bitscores.append(int(info[2].strip()))

    plt.hist(bitscores, 20)
    plt.show()


if __name__ == "__main__":
    plot_bitscores(argv[1])
