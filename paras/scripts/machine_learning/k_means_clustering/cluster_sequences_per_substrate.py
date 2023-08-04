from argparse import ArgumentParser
import os
from random import shuffle

from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from kneed import KneeLocator
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np


from paras.scripts.feature_extraction.sequence_feature_extraction.seq_to_features import get_sequence_features
from paras.scripts.parsers.fasta import read_fasta
from paras.scripts.parsers.parsers import parse_domain_list
from paras.scripts.parsers.iterate_over_dir import iterate_over_dir


def parse_arguments():
    parser = ArgumentParser()
    parser.add_argument('-f', type=str, required=True, help='Input fasta.')
    parser.add_argument('-d', type=str, required=True, help='Domain list directory.')
    parser.add_argument('-o', type=str, required=True, help="Output directory.")

    args = parser.parse_args()
    return args


def run():
    colors = list(cm.tab20(np.linspace(0, 1, 20)))
    shuffle(colors)
    args = parse_arguments()
    if not os.path.exists(args.o):
        os.mkdir(args.o)

    domain_to_seq = read_fasta(args.f)
    for spec, file_path in iterate_over_dir(args.d, '.txt'):
        out_file = os.path.join(args.o, f"{spec}.txt")
        domain_list = parse_domain_list(file_path)
        seq_feature_list = []

        for domain in domain_list:
            sequence = domain_to_seq[domain]
            seq_features = get_sequence_features(sequence)
            seq_feature_list.append(seq_features)

        wccs = []
        nr_clusters = []

        for i in range(1, min(21, len(domain_list) + 1)):
            kmeans = KMeans(n_clusters=i, init='k-means++', n_init=1, random_state=25051989)
            kmeans.fit(seq_feature_list)
            wccs.append(kmeans.inertia_)
            nr_clusters.append(i)

        kneedle = KneeLocator(nr_clusters, wccs, S=1.0, curve="convex", direction="decreasing")
        if kneedle.knee is None:
            nr_clusters = 2
        else:
            nr_clusters = int(kneedle.knee)
        kmeans = KMeans(n_clusters=nr_clusters, init='k-means++', n_init=1, random_state=25051989)
        kmeans.fit(seq_feature_list)

        pca = PCA(2, random_state=25051989)

        vectors = pca.fit_transform(seq_feature_list)

        print(f"Knee for {spec} at:", nr_clusters)
        out_knee = os.path.join(args.o, f"{spec}_elbow.png")
        kneedle.plot_knee()
        plt.savefig(out_knee)
        plt.clf()
        out_png = os.path.join(args.o, f"{spec}.png")

        color_list = []

        with open(out_file, 'w') as out:

            for i, label in enumerate(kmeans.labels_):
                color_list.append(colors[label])
                domain = domain_list[i]
                out.write(f"{domain}\t{label}\n")

        x = []
        y = []

        for vector in vectors:
            x.append(vector[0])
            y.append(vector[1])
        plt.scatter(x, y, color=color_list)
        plt.savefig(out_png)



if __name__ == "__main__":
    run()










