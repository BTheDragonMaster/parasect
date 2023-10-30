import argparse

from paras.scripts.data_analysis.pca.pca import Pca
from paras.scripts.parsers.parsers import parse_domain_list, parse_specificities
from paras.scripts.parsers.fasta import read_fasta
from paras.scripts.feature_extraction.sequence_feature_extraction.seq_to_features import get_sequence_features_bulk


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, required=True, help="Fasta file containing sequences to be contextualised.")
    parser.add_argument('-m', type=str, required=True, help="Path to PCA model.")
    parser.add_argument('-f', type=str, required=True, help="Path to fasta file.")
    parser.add_argument('-d', type=str, default=None, help="Path to domain list of domains that provide context")
    parser.add_argument('-o', type=str, required=True, help="Path to output file")
    parser.add_argument('-s', type=str, required=True, help="Path to specificities file")

    args = parser.parse_args()
    return args


def run():
    args = parse_arguments()

    domain_list = parse_domain_list(args.d)
    domain_to_spec = parse_specificities(args.s)
    domain_to_sequence = read_fasta(args.f)
    for domain in list(domain_to_sequence.keys()):
        if domain not in domain_list:
            del domain_to_sequence[domain]

    domain_to_features, categories = get_sequence_features_bulk(domain_to_sequence)
    domain_to_new = read_fasta(args.i)
    domain_to_new_features, categories = get_sequence_features_bulk(domain_to_new)
    domain_to_features.update(domain_to_new_features)

    vectors = []
    substrate_labels = []
    domain_labels = []

    for domain_id, features in domain_to_features.items():
        vectors.append(features)
        domain_labels.append(domain_id)
        if domain_id in domain_to_spec:
            substrate_labels.append('|'.join(domain_to_spec[domain_id]))
        else:
            substrate_labels.append("Unknown")

    print("Loading model..")

    pca = Pca.from_model(args.m)

    pca.data.set_labels_from_lists(substrate_labels, domain_labels)
    print("Running model..")
    pca.apply(vectors)
    pca.data.write_transformed_vectors(args.o)


if __name__ == "__main__":
    run()
