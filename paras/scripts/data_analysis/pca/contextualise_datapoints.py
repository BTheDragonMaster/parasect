import argparse
import os

import numpy as np

from paras.scripts.data_analysis.pca.pca import Pca
from paras.scripts.parsers.parsers import yield_pocket_features, parse_domain_list, parse_specificities, parse_pocket_features

s = "AAC45929.1.A3|O30408.A3"


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, default=None, help="Path to pocket file containing datapoint(s) to be contextualised.")
    parser.add_argument('-n', type=str, nargs="*", default=None, help="Names of A domains in voxel dataset to contextualise.")
    parser.add_argument('-m', type=str, required=True, help="Path to PCA model.")
    parser.add_argument('-p', type=str, default=None, help="Path to pocket file.")
    parser.add_argument('-d', type=str, default=None, help="Path to domain list of domains that provide context")
    parser.add_argument('-o', type=str, required=True, help="Path to output file")
    parser.add_argument('-s', type=str, required=True, help="Path to specificities file")

    args = parser.parse_args()
    return args


def run():
    args = parse_arguments()

    assert args.i or args.n

    domain_list = parse_domain_list(args.d)
    domain_to_spec = parse_specificities(args.s)
    domain_to_pocket = {}

    for domain in args.n:
        if domain not in domain_list:
            domain_list.append(domain)

    print("Loading data..")
    for domain, pocket, categories in yield_pocket_features(args.p):
        if domain in domain_list:
            domain_to_pocket[domain] = pocket

    for domain, pocket, categories in yield_pocket_features(args.i):
        print("hello", domain)
        domain_to_pocket[domain] = pocket

    vectors = []
    substrate_labels = []
    domain_labels = []

    for domain in domain_to_pocket:
        vectors.append(np.array(domain_to_pocket[domain], dtype='uint8'))
        if domain in domain_to_spec:
            substrate_labels.append('|'.join(domain_to_spec[domain]))
        else:
            substrate_labels.append("Unknown")
        domain_labels.append(domain)

    print("Loading model..")

    pca = Pca.from_model(args.m)

    pca.data.set_labels_from_lists(substrate_labels, domain_labels)
    print("Running model..")
    pca.apply(vectors)
    pca.data.write_transformed_vectors(args.o)


if __name__ == "__main__":
    run()
