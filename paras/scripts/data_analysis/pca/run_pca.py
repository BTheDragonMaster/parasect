import argparse
import os

from sklearn.decomposition import PCA
import numpy as np

from paras.scripts.data_analysis.pca.pca import Pca
from paras.scripts.parsers.parsers import parse_pocket_features, parse_domain_list, parse_specificities
from paras.scripts.parsers.fasta import read_fasta
from paras.scripts.feature_extraction.sequence_feature_extraction.seq_to_features import get_sequence_features_bulk


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', type=str, default=None, help="Path to pocket file.")
    parser.add_argument('-s', type=str, default=None, help="Path to specificities file")
    parser.add_argument('-d', type=str, default=None, help="Path to domain list")
    parser.add_argument('-o', type=str, required=True, help="Path to output directory")
    parser.add_argument('-f', type=str, default=None, help="Path to fasta file")
    parser.add_argument('-m', type=str, default='structure',
                        help="Mode. If mode is 'structure', provide '-p'. If mode is 'sequence', provide '-f'")

    args = parser.parse_args()
    return args


def run():
    args = parse_arguments()

    if not os.path.exists(args.o):
        os.mkdir(args.o)

    print("Loading data..")

    if args.m == 'structure':
        assert args.p
        domain_to_features, categories = parse_pocket_features(args.p, return_categories=True)
    elif args.m =='sequence':
        assert args.f
        domain_to_seq = read_fasta(args.f)
        domain_to_features, categories = get_sequence_features_bulk(domain_to_seq)

    else:
        raise ValueError(f"Expected 'structure' or 'sequence' for mode argument (-m). Got {args.f}.")

    domain_list = parse_domain_list(args.d)
    domain_to_spec = parse_specificities(args.s)

    vectors = []
    substrate_labels = []
    domain_labels = []

    for domain in domain_list:
        vectors.append(np.array(domain_to_features[domain], dtype='uint8'))
        substrate_labels.append('|'.join(domain_to_spec[domain]))
        domain_labels.append(domain)

    print("Processing data..")

    print("Assessing PCA..")
    pca = Pca(PCA(n_components=30))
    pca.data.set_labels_from_lists(substrate_labels, domain_labels)
    pca.fit(vectors)
    knee = pca.find_optimal_nr_components(args.o)

    print("Building PCA..")
    pca = Pca(PCA(n_components=knee))
    pca.data.set_labels_from_lists(substrate_labels, domain_labels)
    pca.fit(vectors)
    pca.save(os.path.join(args.o, "pca.model"))

    pca.apply(vectors)
    pca.data.write_transformed_vectors(os.path.join(args.o, "pcs.txt"))
    pca.write_importance_per_component(categories, args.o)


if __name__ == "__main__":
    run()
