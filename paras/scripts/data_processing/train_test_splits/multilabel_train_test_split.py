#!/usr/bin/env python

import os
import argparse

import numpy as np
from iterstrat.ml_stratifiers import MultilabelStratifiedKFold, MultilabelStratifiedShuffleSplit
from sklearn.preprocessing import MultiLabelBinarizer

from paras.scripts.parsers.parsers import parse_specificities, parse_domain_list
from paras.scripts.data_analysis.count_substrates import count_substrates_from_file, count_substrates
from paras.scripts.data_processing.relabel_data import relabel_from_inclusion_limit


def parse_arguments():
    parser = argparse.ArgumentParser(description="Script to split train, test and validation sets.")
    parser.add_argument('-i', type=str, required=True, help='Path to file containing list of domains.')
    parser.add_argument('-s', type=str, required=True, help='Path to file containing \
    substrate specificities.')
    parser.add_argument('-o', type=str, required=True, help="Output directory.")
    parser.add_argument('-t', type=float, default=0.25, help="Size of test set")
    parser.add_argument('-f', type=int, default=5, help="Fold cross-validation")
    parser.add_argument('-l', type=int, default=11, help="Minimum number of samples required to be considered for training")

    args = parser.parse_args()
    return args


def binarise_data(labels):
    mlb = MultiLabelBinarizer()
    new_labels = mlb.fit_transform(labels)
    return new_labels, mlb.classes_


def assess_label_distribution(all_labels, train_idx, test_idx):
    train_labels = []
    test_labels = []
    for idx in train_idx:
        train_labels += all_labels[idx]
    for idx in test_idx:
        test_labels += all_labels[idx]

    train_counts = count_substrates(train_labels)
    test_counts = count_substrates(test_labels)

    for substrate, train_count in train_counts.items():
        try:
            test_count = test_counts[substrate]
            print(f"{substrate}: {test_count / (train_count + test_count)}")
        except KeyError:
            print(f"{substrate}: 0.0")


def run():
    args = parse_arguments()
    if not os.path.exists(args.o):
        os.mkdir(args.o)

    crossval_path = os.path.join(args.o, 'crossvalidation')
    if not os.path.exists(crossval_path):
        os.mkdir(crossval_path)

    substrate_counts = count_substrates_from_file(args.i, args.s)

    domain_list = parse_domain_list(args.i)
    domain_to_substrates = parse_specificities(args.s)

    domain_to_filtered = relabel_from_inclusion_limit(domain_list, domain_to_substrates, substrate_counts, args.l)

    label_sets = []
    domains = []

    for domain, substrates in domain_to_filtered.items():
        label_sets.append(substrates)
        domains.append(domain)

    binary_labels, labels = binarise_data(label_sets)
    domains = np.array(domains)

    train_test_splitter = MultilabelStratifiedShuffleSplit(n_splits=1, test_size=args.t, random_state=25051989)

    domains_train = []
    domains_test = []

    idx_train = []
    idx_test = []

    labels_train = []
    labels_test = []

    binary_train_labels = []

    for train_x, test_x in train_test_splitter.split(domains, binary_labels):

        for idx in train_x:
            domains_train.append(domains[idx])
            idx_train.append(idx)
            labels_train.append(label_sets[idx])
            binary_train_labels.append(binary_labels[idx])

        for idx in test_x:
            domains_test.append(domains[idx])
            idx_test.append(idx)
            labels_test.append(label_sets[idx])

        break

    assess_label_distribution(label_sets, idx_train, idx_test)

    crossval_sets = make_crossval_sets(domains_train, binary_train_labels, args.f)

    for i, (crossval_train, crossval_test) in enumerate(crossval_sets):
        print(f"Crossval set: {i}")
        assess_label_distribution(labels_train, crossval_train, crossval_test)
        print("\n")

        file_path_train = os.path.join(crossval_path, f'crossval_train_{i:02d}.txt')
        file_path_test = os.path.join(crossval_path, f'crossval_test_{i:02d}.txt')

        with open(file_path_train, 'w') as train_file:
            for index in crossval_train:
                domain = domains_train[index]
                train_file.write(f"{domain}\n")

        with open(file_path_test, 'w') as test_file:
            for index in crossval_test:
                domain = domains_train[index]
                test_file.write(f"{domain}\n")

    train_file = os.path.join(args.o, "train.txt")
    file_path_included = os.path.join(args.o, f'included_substrates.txt')

    with open(file_path_included, 'w') as included:
        labels.sort()

        for label in labels:
            included.write(f"{label}\n")

    with open(train_file, 'w') as train:
        for domain in domains_train:
            train.write(f"{domain}\n")

    test_file = os.path.join(args.o, "test.txt")
    with open(test_file, 'w') as test:
        for domain in domains_test:
            test.write(f"{domain}\n")


def make_crossval_sets(train_x, train_y, fold_validation):
    crossval_sets = []
    stratifier = MultilabelStratifiedKFold(n_splits=fold_validation, shuffle=True, random_state=25051989)
    for train_idx, test_idx in stratifier.split(train_x, train_y):
        crossval_sets.append((train_idx, test_idx))
    return crossval_sets


if __name__ == "__main__":
    run()
