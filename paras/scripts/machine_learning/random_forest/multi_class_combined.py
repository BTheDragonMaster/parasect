from paras.scripts.parsers.parsers import parse_specificities, parse_domain_list, parse_pocket_features, \
    parse_stach_codes
from paras.scripts.parsers.tabular import Tabular
from paras.scripts.data_processing.get_domain_list import make_domain_mapping
from paras.scripts.feature_extraction.sequence_feature_extraction.seq_to_features import get_sequence_features, PROPERTIES_FILE
from paras.scripts.machine_learning.random_forest.multi_class_no_moments import plot_amino_acid_heatmap
from joblib import dump, load

from sklearn.ensemble import RandomForestClassifier
import argparse
from typing import Optional
from collections import OrderedDict
import os


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', type=str, required=True, help="Path to specificities file")
    parser.add_argument('-d', type=str, required=True, help="Path to domain file")
    parser.add_argument('-p', type=str, required=True, help="Path to pocket file.")
    parser.add_argument('-o', type=str, default=None, help="Path to output directory")
    parser.add_argument('-c', type=str, default=None, help="Path to classifier directory")
    parser.add_argument('-fa', type=str, default=None, help="Path to input tabular containing constant-length specificity-conferring codes")
    parser.add_argument('-code', type=str, default='34', help="Sequence mode. 'stach' or '34'")
    parser.add_argument('-n', type=int, default=1000, help="Nr of trees in RF")
    parser.add_argument('-t', type=int, default=None, help="Nr of threads. -1 uses all available cores")
    parser.add_argument('-train', type=str, default=None, help="Path to train domain file")
    parser.add_argument('-test', type=str, default=None, help="Path to test domain file")
    parser.add_argument('-sampling', type=str, default=None, help="balanced, balanced_subsample, or None")

    args = parser.parse_args()
    return args


class RandomForestParasect:
    def __init__(self, pocket_file, sequence_file, specificities_file, domain_file, sequence_mode='34'):
        print("Parsing data..")
        sequence_feature_names = Tabular(PROPERTIES_FILE, [0]).categories[1:]
        domain_to_stach, domain_to_34 = parse_stach_codes(sequence_file)

        self.sequence_features = []

        if sequence_mode == '34':
            self.domain_to_sequence = domain_to_34
            nr_res = 34
        elif sequence_mode == 'stach':
            self.domain_to_sequence = domain_to_stach
            nr_res = 9
        else:
            raise ValueError("Sequence mode must be either 'stach' or '34'.")

        for i in range(nr_res):
            for sequence_feature_name in sequence_feature_names:
                self.sequence_features.append(f"res{i}_{sequence_feature_name}")

        domain_mapping = make_domain_mapping(domain_file)
        old_domain_to_pocket_vector = parse_pocket_features(pocket_file)
        self.domain_to_pocket_vector = {}
        for old_enzyme, vector in old_domain_to_pocket_vector.items():
            enzyme = domain_mapping[old_enzyme]
            self.domain_to_pocket_vector[enzyme] = vector

        self.domain_to_substrates = parse_specificities(specificities_file)
        with open(pocket_file, 'r') as pocket_info:
            self.pocket_features = pocket_info.readline().split('\t')[1:]

        self.features = self.pocket_features + self.sequence_features

    def train_model(self, train_domain_file, out_folder, sampling_method, n_jobs: Optional[int] = None,
                    n_estimators: int = 1000, oob_score: bool = True):
        print("Building data..")

        train_x, train_y, substrate_order, domain_order = self.build_data(train_domain_file)

        seq_train_x, seq_train_y, substrate_order, domain_order = self.build_data(train_domain_file, mode='sequence')

        vox_train_x, vox_train_y, substrate_order, domain_order = self.build_data(train_domain_file, mode='structure')

        if sampling_method == 'balanced':

            sequence_classifier = RandomForestClassifier(n_estimators=n_estimators, n_jobs=n_jobs, oob_score=oob_score,
                                                random_state=25051989, class_weight='balanced')
            structure_classifier = RandomForestClassifier(n_estimators=n_estimators, n_jobs=n_jobs, oob_score=oob_score,
                                                random_state=25051989, class_weight='balanced')
            combined_classifier = RandomForestClassifier(n_estimators=n_estimators, n_jobs=n_jobs, oob_score=oob_score,
                                                          random_state=25051989, class_weight='balanced')
        elif sampling_method == 'balanced_subsample':
            sequence_classifier = RandomForestClassifier(n_estimators=n_estimators, n_jobs=n_jobs, oob_score=oob_score,
                                                random_state=25051989, class_weight='balanced_subsample')
            structure_classifier = RandomForestClassifier(n_estimators=n_estimators, n_jobs=n_jobs, oob_score=oob_score,
                                                         random_state=25051989, class_weight='balanced_subsample')
            combined_classifier = RandomForestClassifier(n_estimators=n_estimators, n_jobs=n_jobs, oob_score=oob_score,
                                                          random_state=25051989, class_weight='balanced_subsample')
        else:
            sequence_classifier = RandomForestClassifier(n_estimators=n_estimators, n_jobs=n_jobs, oob_score=oob_score,
                                                random_state=25051989)
            structure_classifier = RandomForestClassifier(n_estimators=n_estimators, n_jobs=n_jobs, oob_score=oob_score,
                                                         random_state=25051989)
            combined_classifier = RandomForestClassifier(n_estimators=n_estimators, n_jobs=n_jobs, oob_score=oob_score,
                                                          random_state=25051989)
        print("Training combined classifier..")
        combined_classifier.fit(train_x, train_y)
        print("Training sequence classifier..")
        sequence_classifier.fit(seq_train_x, seq_train_y)
        print("Training structure classifier..")
        structure_classifier.fit(vox_train_x, vox_train_y)
        print("Oob scores:")
        print("Sequence:", sequence_classifier.oob_score_)
        print("Structure:", structure_classifier.oob_score_)
        print("Combined:", combined_classifier.oob_score_)

        sequence_out = os.path.join(out_folder, 'sequence_only.parasect')
        structure_out = os.path.join(out_folder, 'structure_only.parasect')
        combined_out = os.path.join(out_folder, 'combined.parasect')
        dump(sequence_classifier, sequence_out)
        dump(structure_classifier, structure_out)
        dump(combined_classifier, combined_out)

        sequence_importance_out = os.path.join(out_folder, 'sequence_only_importances.txt')
        structure_importance_out = os.path.join(out_folder, 'structure_only_importances.txt')
        combined_importance_out = os.path.join(out_folder, 'combined_importances.txt')

        with open(sequence_importance_out, 'w') as feature_file:
            for i, importance in enumerate(sequence_classifier.feature_importances_):
                feature_name = self.sequence_features[i]
                feature_file.write(f"{feature_name}\t{importance}\n")

        with open(structure_importance_out, 'w') as feature_file:
            for i, importance in enumerate(structure_classifier.feature_importances_):
                feature_name = self.pocket_features[i]
                feature_file.write(f"{feature_name}\t{importance}\n")

        with open(combined_importance_out, 'w') as feature_file:
            for i, importance in enumerate(combined_classifier.feature_importances_):
                feature_name = self.features[i]
                feature_file.write(f"{feature_name}\t{importance}\n")


    def test_model(self, test_domain_file, classifier_dir):
        for classifier in os.listdir(classifier_dir):
            if classifier.endswith('.parasect'):
                classifier_file = os.path.join(classifier_dir, classifier)
                classifier_name = classifier.split('.parasect')[0]

                if classifier_name == "sequence_only":

                    test_x, test_y, substrate_order, domain_order = self.build_data(test_domain_file, mode='sequence')
                elif classifier_name == 'structure_only':
                    test_x, test_y, substrate_order, domain_order = self.build_data(test_domain_file, mode='structure')
                elif classifier_name == 'combined':
                    test_x, test_y, substrate_order, domain_order = self.build_data(test_domain_file, mode='combined')
                else:
                    raise ValueError("Name of classifier should be either structure_only, sequence_only or combined")

                classifier = load(classifier_file)
                classes = classifier.classes_
                predict_y = classifier.predict_proba(test_x)

                correct = 0
                incorrect = 0

                substrate_to_correct = {}
                substrate_to_actual = OrderedDict()
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

                print(f"Predictive power for {classifier_name}:")
                print(f"Total points: {len(test_y)}")
                print(f"Accuracy: {correct / (correct + incorrect)}")

                for substrate, correct_list in substrate_to_correct.items():
                    substrate_to_correct[substrate] = correct_list.count(1) / len(correct_list)
                plot_out = os.path.join(classifier_dir, f"{classifier_name}.png")

                plot_amino_acid_heatmap(substrate_to_actual, plot_out)

    def build_data(self, domain_file: str, mode='combined'):
        domains = parse_domain_list(domain_file)

        x = []
        y = []

        substrate_order = []
        domain_order = []

        for domain in domains:
            substrate = self.domain_to_substrates[domain][0]

            sequence_vector = get_sequence_features(self.domain_to_sequence[domain])
            if mode == 'combined':
                pocket_vector = self.domain_to_pocket_vector[domain]
                vector = pocket_vector + sequence_vector
            elif mode == 'sequence':
                vector = sequence_vector
            elif mode == 'structure':
                pocket_vector = self.domain_to_pocket_vector[domain]
                vector = pocket_vector
            else:
                raise ValueError("Data building mode must be either combined, sequence, or structure")
            x.append(vector)
            y.append(substrate)
            domain_order.append(domain)
            substrate_order.append(substrate)

        return x, y, substrate_order, domain_order


def run():
    args = parse_arguments()

    random_forest = RandomForestParasect(args.p, args.fa, args.s, args.d, args.code)
    if args.c:
        assert args.test
        random_forest.test_model(args.test, args.c)
    elif args.o:
        if not os.path.exists(args.o):
            os.mkdir(args.o)
        assert args.train
        random_forest.train_model(args.train, args.o, args.sampling, n_estimators=args.n, n_jobs=args.t, )
        if args.test:
            random_forest.test_model(args.test, args.o)


if __name__ == "__main__":
    run()
