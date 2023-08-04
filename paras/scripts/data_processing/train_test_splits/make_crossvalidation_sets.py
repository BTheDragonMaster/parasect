import argparse
import os

from iterstrat.ml_stratifiers import MultilabelStratifiedKFold
import numpy as np

from paras.scripts.parsers.parsers import parse_specificities, parse_domain_list
from paras.scripts.data_processing.relabel_data import relabel_from_substrate_list
from paras.scripts.data_processing.train_test_splits.multilabel_train_test_split import binarise_data, \
    assess_label_distribution


def parse_arguments():
    parser = argparse.ArgumentParser(description="Script to split train, test and validation sets.")
    parser.add_argument('-d', type=str, required=True, help='Path to file containing list of domains.')
    parser.add_argument('-s', type=str, required=True, help='Path to file containing \
    substrate specificities.')
    parser.add_argument('-i', type=str, required=True, help='Path to file containing list of included substrates.')
    parser.add_argument('-o', type=str, required=True, help="Output directory.")
    parser.add_argument('-f', type=int, default=5, help="Fold cross-validation")

    args = parser.parse_args()
    return args


def run():
    args = parse_arguments()
    if not os.path.exists(args.o):
        os.mkdir(args.o)

    crossval_path = os.path.join(args.o, 'crossvalidation')
    if not os.path.exists(crossval_path):
        os.mkdir(crossval_path)

    domain_list = parse_domain_list(args.d)
    domain_to_substrates = parse_specificities(args.s)

    domain_to_filtered = relabel_from_substrate_list(domain_list, domain_to_substrates, args.i)

    label_sets = []
    domains = []

    for domain, substrates in domain_to_filtered.items():
        label_sets.append(substrates)
        domains.append(domain)

    binary_labels, labels = binarise_data(label_sets)
    domains = np.array(domains)

    crossval_sets = make_crossval_sets(domains, binary_labels, args.f)

    for i, (crossval_train, crossval_test) in enumerate(crossval_sets):
        print(f"Crossval set: {i}")
        assess_label_distribution(label_sets, crossval_train, crossval_test)
        print("\n")

        file_path_train = os.path.join(crossval_path, f'crossval_train_{i:02d}.txt')
        file_path_test = os.path.join(crossval_path, f'crossval_test_{i:02d}.txt')

        with open(file_path_train, 'w') as train_file:
            for index in crossval_train:
                domain = domains[index]
                train_file.write(f"{domain}\n")

        with open(file_path_test, 'w') as test_file:
            for index in crossval_test:
                domain = domains[index]
                test_file.write(f"{domain}\n")


def make_crossval_sets(train_x, train_y, fold_validation):
    crossval_sets = []
    stratifier = MultilabelStratifiedKFold(n_splits=fold_validation, shuffle=True, random_state=25051989)
    for train_idx, test_idx in stratifier.split(train_x, train_y):
        crossval_sets.append((train_idx, test_idx))
    return crossval_sets


if __name__ == "__main__":
    run()