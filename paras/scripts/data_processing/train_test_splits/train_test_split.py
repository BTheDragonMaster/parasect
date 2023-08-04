#!/usr/bin/env python

import os
import argparse

from paras.scripts.parsers.parsers import parse_specificities, parse_domain_list
from sklearn.model_selection import train_test_split, StratifiedKFold


def parse_arguments():
    parser = argparse.ArgumentParser(description="Script to split train, test and validation sets.")
    parser.add_argument('-i', '--domains', type=str, required=True, help='Path to file containing list of domains.')
    parser.add_argument('-s', '--substrates', type=str, required=True, help='Path to file containing \
    substrate specificities.')
    parser.add_argument('-o', type=str, required=True, help="Output directory.")
    parser.add_argument('-t', type=float, default=0.1, help="Size of test set")
    parser.add_argument('-f', type=int, default=10, help="Fold cross-validation")
    parser.add_argument('-l', type=int, default=20, help="Minimum number of samples required to be considered for training")

    args = parser.parse_args()
    return args


def run():
    args = parse_arguments()
    if not os.path.exists(args.o):
        os.mkdir(args.o)

    crossval_path = os.path.join(args.o, 'crossvalidation')
    if not os.path.exists(crossval_path):
        os.mkdir(crossval_path)

    domain_list = parse_domain_list(args.domains)
    substrates = []
    all_substrates = []
    domain_to_substrates = parse_specificities(args.substrates)

    for domain in domain_list:
        substrate = domain_to_substrates[domain][0]
        substrates.append(substrate)
        all_substrates += domain_to_substrates[domain]

    substrate_counts = {}

    for substrate in all_substrates:
        if substrate not in substrate_counts:
            substrate_counts[substrate] = 0
        substrate_counts[substrate] += 1

    limit = args.l

    filtered_domains = []
    filtered_substrates = []

    remaining_domains = []
    remaining_substrates = []

    for i, substrate in enumerate(substrates):
        if substrate_counts[substrate] >= limit:
            filtered_substrates.append(substrate)
            filtered_domains.append(domain_list[i])
        else:
            remaining_domains.append(domain_list[i])
            remaining_substrates.append(substrate)

    print(len(filtered_domains))
    print(len(filtered_substrates))

    train_x, test_x, train_y, test_y = train_test_split(filtered_domains, filtered_substrates, test_size=args.t,
                                                        stratify=filtered_substrates, random_state=25051989)

    crossval_sets = make_crossval_sets(train_x, train_y, args.f)

    for i, (crossval_train, crossval_test) in enumerate(crossval_sets):
        file_path_train = os.path.join(crossval_path, f'crossval_train_{i:02d}.txt')
        file_path_test = os.path.join(crossval_path, f'crossval_test_{i:02d}.txt')
        with open(file_path_train, 'w') as train_file:
            for index in crossval_train:
                domain = filtered_domains[index]
                train_file.write(f"{domain}\n")

        with open(file_path_test, 'w') as test_file:
            for index in crossval_test:
                domain = filtered_domains[index]
                test_file.write(f"{domain}\n")

    train_file = os.path.join(args.o, "train.txt")
    with open(train_file, 'w') as train:
        for domain in train_x:
            train.write(f"{domain}\n")

    test_file = os.path.join(args.o, "test.txt")
    with open(test_file, 'w') as test:
        for domain in test_x:
            test.write(f"{domain}\n")

    filtered_file = os.path.join(args.o, "filtered.txt")
    with open(filtered_file, 'w') as filtered:
        for domain in remaining_domains:
            filtered.write(f"{domain}\n")


def make_crossval_sets(train_x, train_y, fold_validation):
    stratifier = StratifiedKFold(n_splits=fold_validation, random_state=25051989, shuffle=True)
    splits = stratifier.split(train_x, train_y)
    return splits


if __name__ == "__main__":
    run()
