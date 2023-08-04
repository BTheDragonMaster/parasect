import argparse
import os
from paras.scripts.parsers.parsers import parse_interactions, parse_domain_list
from sklearn.model_selection import train_test_split, StratifiedKFold


def parse_arguments():
    parser = argparse.ArgumentParser(description="Script to split train, test and validation sets per substrate.")
    parser.add_argument('-i', type=str, required=True, help='Path to folder containing files storing interaction per substrate.')
    parser.add_argument('-o', type=str, required=True, help="Output directory.")
    parser.add_argument('-d', type=str, required=True, help="Path to domain list.")
    parser.add_argument('-t', type=float, default=0.1, help="Test size")
    parser.add_argument('-f', type=int, default=10, help="Fold crossvalidation")
    parser.add_argument('-l', type=int, default=20, help="Minimum number of positive examples for a substrate to build a predictor.")

    args = parser.parse_args()
    return args


def run():
    args = parse_arguments()
    if not os.path.exists(args.o):
        os.mkdir(args.o)

    counter = 0
    untrained_counter = 0
    domain_list = set(parse_domain_list(args.d))

    for substrate_file in os.listdir(args.i):
        domains = []
        interactions = []
        substrate_path = os.path.join(args.i, substrate_file)
        if substrate_path.endswith('.txt') and os.path.isfile(substrate_path):
            substrate_name = substrate_file.split('.txt')[0]
            substrate_out = os.path.join(args.o, substrate_name)

            domain_to_interaction = parse_interactions(substrate_path)
            for domain, interaction in domain_to_interaction.items():
                if domain in domain_list:
                    domains.append(domain)
                    interactions.append(interaction)

            true_counts = interactions.count(1)
            if true_counts >= args.l:

                if not os.path.exists(substrate_out):
                    os.mkdir(substrate_out)

                crossval_path = os.path.join(substrate_out, 'crossvalidation')
                if not os.path.exists(crossval_path):
                    os.mkdir(crossval_path)

                counter += 1
                print(f"{substrate_name} has enough examples with {true_counts} data points.")
                train_x, test_x, train_y, test_y = train_test_split(domains, interactions, test_size=args.t,
                                                                    stratify=interactions, random_state=25051989)

                crossval_sets = make_crossval_sets(train_x, train_y, args.f)

                for i, (crossval_train, crossval_test) in enumerate(crossval_sets):
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

                train_file = os.path.join(substrate_out, "train.txt")
                with open(train_file, 'w') as train:
                    for domain in train_x:
                        train.write(f"{domain}\n")

                test_file = os.path.join(substrate_out, "test.txt")
                with open(test_file, 'w') as test:
                    for domain in test_x:
                        test.write(f"{domain}\n")

            else:
                untrained_counter += 1

    print(f"Training sets made for {counter} substrates.")
    print(f"No training sets made for {untrained_counter} substrates.")


def make_crossval_sets(train_x, train_y, fold_validation):
    stratifier = StratifiedKFold(n_splits=fold_validation, random_state=25051989, shuffle=True)
    splits = stratifier.split(train_x, train_y)
    return splits


if __name__ == "__main__":
    run()
