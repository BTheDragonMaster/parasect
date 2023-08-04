from argparse import ArgumentParser
import os
from random import shuffle

from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from kneed import KneeLocator
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np

from paras.scripts.parsers.parsers import parse_domain_list, parse_bitscores
from paras.scripts.parsers.iterate_over_dir import iterate_over_dir


def parse_arguments():
    parser = ArgumentParser()
    parser.add_argument('-b', type=str, required=True, help='Bitscore path.')
    parser.add_argument('-d', type=str, required=True, help='Directory of domain lists per substrate.')
    parser.add_argument('-o', type=str, required=True, help="Output directory.")
    parser.add_argument('-a', type=str, required=True, help='Domain list of all domains.')

    args = parser.parse_args()
    return args


def run():
    colors = list(cm.tab20(np.linspace(0, 1, 20)))
    shuffle(colors)
    args = parse_arguments()
    if not os.path.exists(args.o):
        os.mkdir(args.o)

    out_file = os.path.join(args.o, f"all.txt")

    domain_list = parse_domain_list(args.a)
    print("Processing data..")
    domain_to_bitscore_vector = parse_bitscores(args.b, domain_list)
    print("Processed data")
    vectors = []

    for domain in domain_list:
        vectors.append(domain_to_bitscore_vector[domain])

    wccs = []
    nr_clusters = []

    for i in range(1, min(21, len(domain_list) + 1)):
        print(f"Kmeans iteration {i} for all")
        kmeans = KMeans(n_clusters=i, init='k-means++', n_init=1, random_state=25051989)
        kmeans.fit(vectors)
        wccs.append(kmeans.inertia_)
        nr_clusters.append(i)

    kneedle = KneeLocator(nr_clusters, wccs, S=1.0, curve="convex", direction="decreasing")
    if kneedle.knee is None:
        nr_clusters = 2
    else:
        nr_clusters = int(kneedle.knee)
    kmeans = KMeans(n_clusters=nr_clusters, init='k-means++', n_init=1, random_state=25051989)
    kmeans.fit(vectors)

    pca = PCA(2, random_state=25051989)

    vectors = pca.fit_transform(vectors)

    print(f"Nr of clusters for all:", nr_clusters)
    out_knee = os.path.join(args.o, f"all_elbow.png")
    kneedle.plot_knee()
    plt.savefig(out_knee)
    plt.clf()
    plt.close()
    out_png = os.path.join(args.o, f"all.png")

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
    plt.clf()
    plt.close()

    for spec, file_path in iterate_over_dir(args.d, '.txt'):
        out_file = os.path.join(args.o, f"{spec}.txt")
        domain_list = parse_domain_list(file_path)

        vectors = []

        for domain in domain_list:
            vectors.append(domain_to_bitscore_vector[domain])

        wccs = []
        nr_clusters = []

        for i in range(1, min(21, len(domain_list) + 1)):
            print(f"Kmeans iteration {i} for {spec}")
            kmeans = KMeans(n_clusters=i, init='k-means++', n_init=1, random_state=25051989)
            kmeans.fit(vectors)
            wccs.append(kmeans.inertia_)
            nr_clusters.append(i)

        kneedle = KneeLocator(nr_clusters, wccs, S=1.0, curve="convex", direction="decreasing")
        if kneedle.knee is None:
            nr_clusters = 2
        else:
            nr_clusters = int(kneedle.knee)
        kmeans = KMeans(n_clusters=nr_clusters, init='k-means++', n_init=1, random_state=25051989)
        kmeans.fit(vectors)

        pca = PCA(2, random_state=25051989)

        vectors = pca.fit_transform(vectors)

        print(f"Nr of clusters for {spec}:", nr_clusters)
        out_knee = os.path.join(args.o, f"{spec}_elbow.png")
        kneedle.plot_knee()
        plt.savefig(out_knee)
        plt.clf()
        plt.close()
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
        plt.clf()
        plt.close()




if __name__ == "__main__":
    run()

