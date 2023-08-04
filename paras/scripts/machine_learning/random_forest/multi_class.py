from paras.scripts.parsers.parsers import parse_moment_vectors, parse_specificities, parse_domain_list
from paras.scripts.data_processing.get_domain_list import make_domain_mapping
from joblib import dump, load

from sklearn.ensemble import RandomForestClassifier
import argparse
from typing import Optional
from pprint import pprint
from collections import OrderedDict


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', type=str, required=True, help="Path to specificities file")
    parser.add_argument('-m', type=str, required=True, help="Path to moment file.")
    parser.add_argument('-d', type=str, required=True, help="Path to domain file")
    parser.add_argument('-o', type=str, default=None, help="Path to output classifier")
    parser.add_argument('-c', type=str, default=None, help="Path to classifier")
    parser.add_argument('-n', type=int, default=1000, help="Nr of trees in RF")
    parser.add_argument('-t', type=int, default=None, help="Nr of threads. -1 uses all available cores")
    parser.add_argument('-train', type=str, default=None, help="Path to train domain file")
    parser.add_argument('-test', type=str, default=None, help="Path to test domain file")
    parser.add_argument('-sampling', type=str, default=None, help="balanced, balanced_subsample, or None")

    args = parser.parse_args()
    return args


class RandomForestParasect:
    def __init__(self, moment_file, specificities_file, domain_file):
        domain_mapping = make_domain_mapping(domain_file)
        old_domain_to_moment_vector = parse_moment_vectors(moment_file)
        self.domain_to_moment_vector = {}
        for old_enzyme, vector in old_domain_to_moment_vector.items():
            enzyme = domain_mapping[old_enzyme]
            self.domain_to_moment_vector[enzyme] = vector

        self.domain_to_substrates = parse_specificities(specificities_file)

    def train_model(self, train_domain_file, out_path, sampling_method, n_jobs: Optional[int] = None,
                    n_estimators: int = 1000, oob_score: bool = True):
        print("Building data..")

        train_x, train_y, substrate_order, domain_order = self.build_data(train_domain_file)

        if sampling_method == 'balanced':

            classifier = RandomForestClassifier(n_estimators=n_estimators, n_jobs=n_jobs, oob_score=oob_score,
                                                random_state=25051989, class_weight='balanced')
        elif sampling_method == 'balanced_subsample':
            classifier = RandomForestClassifier(n_estimators=n_estimators, n_jobs=n_jobs, oob_score=oob_score,
                                                random_state=25051989, class_weight='balanced_subsample')
        else:
            classifier = RandomForestClassifier(n_estimators=n_estimators, n_jobs=n_jobs, oob_score=oob_score,
                                                random_state=25051989)
        print("Training classifier..")
        classifier.fit(train_x, train_y)
        print(classifier.oob_score_)
        dump(classifier, out_path)

    def test_model(self, test_domain_file, classifier_file):
        test_x, test_y, substrate_order, domain_order = self.build_data(test_domain_file)
        classifier = load(classifier_file)
        classes = classifier.classes_
        predict_y = classifier.predict_proba(test_x)

        correct = 0
        incorrect = 0

        substrate_to_correct = {}
        substrate_to_actual = {}
        substrate_to_rates = {}

        for i, true_y in enumerate(test_y):
            probs_and_substrates = []
            probabilities = predict_y[i]
            for j, substrate in enumerate(classes):
                probability = probabilities[j]
                probs_and_substrates.append((probability, substrate))

            probs_and_substrates.sort(reverse=True)

            best_substrate = probs_and_substrates[0][1]
            if true_y not in substrate_to_correct:
                substrate_to_correct[true_y] = []
                substrate_to_actual[true_y] = {}
            if true_y not in substrate_to_rates:
                substrate_to_rates[true_y] = OrderedDict({"TP": 0,
                                                          "FP": 0,
                                                          "FN": 0})
            if best_substrate not in substrate_to_rates:
                substrate_to_rates[best_substrate] = OrderedDict({"TP": 0,
                                                                  "FP": 0,
                                                                  "FN": 0})

            if best_substrate not in substrate_to_actual[true_y]:
                substrate_to_actual[true_y][best_substrate] = 0
            substrate_to_actual[true_y][best_substrate] += 1

            if best_substrate == true_y:
                correct += 1
                substrate_to_correct[true_y].append(1)
                substrate_to_rates[true_y]["TP"] += 1

            else:
                incorrect += 1
                substrate_to_correct[true_y].append(0)
                substrate_to_rates[best_substrate]["FP"] += 1
                substrate_to_rates[true_y]["FN"] += 1

        print(f"Total points: {len(test_y)}")
        print(f"Accuracy: {correct / (correct + incorrect)}")

        for substrate, correct_list in substrate_to_correct.items():
            substrate_to_correct[substrate] = correct_list.count(1) / len(correct_list)

        pprint(substrate_to_correct)
        pprint(substrate_to_actual)
        pprint(substrate_to_rates)

    def build_data(self, domain_file: str, under_sample=False, over_sample=False):
        domains = parse_domain_list(domain_file)
        if under_sample:
            assert not over_sample
        if over_sample:
            assert not under_sample

        x = []
        y = []

        substrate_order = []
        domain_order = []

        for domain in domains:
            substrate = self.domain_to_substrates[domain][0]
            x.append(self.domain_to_moment_vector[domain])
            y.append(substrate)
            domain_order.append(domain)
            substrate_order.append(substrate)

        return x, y, substrate_order, domain_order


def run():
    args = parse_arguments()

    random_forest = RandomForestParasect(args.m, args.s, args.d)
    if args.c:
        assert args.test
        random_forest.test_model(args.test, args.c)
    elif args.o:
        assert args.train
        random_forest.train_model(args.train, args.o, n_estimators=args.n, n_jobs=args.t, sampling_method=args.sampling)
        if args.test:
            random_forest.test_model(args.test, args.o)


if __name__ == "__main__":
    run()
