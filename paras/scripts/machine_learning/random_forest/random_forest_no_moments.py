from paras.scripts.parsers.parsers import parse_pocket_features, parse_bitvectors, parse_specificities, \
    parse_domain_list
from paras.scripts.data_processing.get_domain_list import make_domain_mapping
from joblib import dump, load

from sklearn.ensemble import RandomForestClassifier
from imblearn.under_sampling import RandomUnderSampler
from imblearn.over_sampling import SMOTE
import argparse
from statistics import mean, median
from typing import Optional
import numpy as np


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', type=str, required=True, help="Path to bitvector file.")
    parser.add_argument('-p', type=str, required=True, help="Path to pocket file.")
    parser.add_argument('-s', type=str, required=True, help="Path to specificities file")
    parser.add_argument('-d', type=str, required=True, help="Path to domain file")
    parser.add_argument('-o', type=str, default=None, help="Path to output classifier")
    parser.add_argument('-c', type=str, default=None, help="Path to classifier")
    parser.add_argument('-n', type=int, default=100, help="Nr of trees in RF")
    parser.add_argument('-t', type=int, default=None, help="Nr of threads. -1 uses all available cores")
    parser.add_argument('-train', type=str, default=None, help="Path to train domain file")
    parser.add_argument('-test', type=str, default=None, help="Path to test domain file")
    parser.add_argument('-sampling', type=str, default=None, help="balanced, balanced_subsample, over_sample, under_sample or None")

    args = parser.parse_args()
    return args


class RFTrainingPoint:
    def __init__(self, domain: str, substrate: str, interaction: int):
        self.domain = domain
        self.substrate = substrate
        self.interaction = interaction
        self.vector = []

    def build_vector(self, pocket_mapping, bitvector_mapping):
        moment_vector = pocket_mapping[self.domain]
        bitvector = bitvector_mapping[self.substrate]
        self.vector = moment_vector + bitvector
        if 'branched' in self.substrate:
            self.vector.append(1)
        else:
            self.vector.append(0)


class RandomForestParasect:
    def __init__(self, bitvector_file, pocket_file, specificities_file, domain_file):
        domain_mapping = make_domain_mapping(domain_file)
        print("Parsing data...")
        old_domain_to_pocket_vector = parse_pocket_features(pocket_file)
        self.domain_to_pocket_vector = {}
        for old_enzyme, vector in old_domain_to_pocket_vector.items():
            enzyme = domain_mapping[old_enzyme]
            self.domain_to_pocket_vector[enzyme] = vector

        self.substrate_to_bitvector = parse_bitvectors(bitvector_file)
        self.domain_to_substrates = parse_specificities(specificities_file)

    def train_model(self, train_domain_file, out_path, sampling_method, n_jobs: Optional[int] = None,
                    n_estimators: int = 1000, oob_score: bool = True):
        print("Building data..")
        if sampling_method == 'over_sample':
            train_x, train_y, substrate_order, domain_order = self.build_data(train_domain_file, over_sample=True)
        elif sampling_method == 'under_sample':
            train_x, train_y, substrate_order, domain_order = self.build_data(train_domain_file, under_sample=True)
        else:
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

        FP = 0
        FN = 0
        TP = 0
        TN = 0

        probs_fp = []
        probs_fn = []
        probs_tp = []
        probs_tn = []

        fp_to_true = {}
        domain_to_predictions = {}
        domain_to_correct = {}
        domain_to_best_prob = {}

        if classes[0] == 1:
            true_index = 0
            false_index = 1
        else:
            true_index = 1
            false_index = 0

        for i, true_y in enumerate(test_y):
            prob_true = predict_y[i][true_index]
            prob_false = predict_y[i][false_index]
            if prob_true > prob_false:
                predicted_y = 1
            else:
                predicted_y = 0

            tested_substrate = substrate_order[i]
            tested_domain = domain_order[i]
            true_substrates = self.domain_to_substrates[tested_domain]
            if tested_domain not in domain_to_correct:
                domain_to_correct[tested_domain] = False
                domain_to_predictions[tested_domain] = []
                domain_to_best_prob[tested_domain] = 0.0

            if prob_true > domain_to_best_prob[tested_domain]:
                domain_to_best_prob[tested_domain] = prob_true
                if tested_substrate in true_substrates:
                    domain_to_correct[tested_domain] = True
                else:
                    domain_to_correct[tested_domain] = False

            domain_to_predictions[tested_domain].append((prob_true, tested_substrate,
                                                         tested_substrate in true_substrates))

            if true_y:
                if not predicted_y:
                    FN += 1
                    probs_fn.append(prob_true)
                else:
                    TP += 1
                    probs_tp.append(prob_true)
            elif not true_y:
                if not predicted_y:
                    TN += 1
                    probs_tn.append(prob_true)
                else:
                    FP += 1
                    probs_fp.append(prob_true)
                    if tested_substrate not in fp_to_true:
                        fp_to_true[tested_substrate] = {}
                    for true_substrate in true_substrates:
                        if true_substrate not in fp_to_true[tested_substrate]:
                            fp_to_true[tested_substrate][true_substrate] = 0
                        fp_to_true[tested_substrate][true_substrate] += 1

        print(f"Total points: {len(test_y)}")

        print(f"True positive: {TP} (average probability of interaction: {mean(probs_tp)}; median: {median(probs_tp)})")
        print(f"False positive: {FP} (average probability of interaction: {mean(probs_fp)}; median: {median(probs_fp)})")
        print(f"True negative: {TN} (average probability of interaction: {mean(probs_tn)}; median: {median(probs_tn)})")
        print(f"False negative: {FN} (average probability of interaction: {mean(probs_fn)}; median: {median(probs_fn)})")
        # pprint(fp_to_true)
        print(f"Correct domain-level predictions: {list(domain_to_correct.values()).count(True)}")
        print(f"Incorrect domain-level predictions: {list(domain_to_correct.values()).count(False)}")
        print(f"Correct percentage: {list(domain_to_correct.values()).count(True) / (list(domain_to_correct.values()).count(False) + list(domain_to_correct.values()).count(True))} ")
        for domain in domain_to_predictions:
            domain_to_predictions[domain].sort(reverse=True)
        # pprint(domain_to_predictions)

    def build_data(self, domain_file: str, under_sample=False, over_sample=False):
        domains = parse_domain_list(domain_file)
        if under_sample:
            assert not over_sample
        if over_sample:
            assert not under_sample

        bitvector_length = len(self.substrate_to_bitvector[list(self.substrate_to_bitvector.keys())[0]])
        pocket_vector_length = len(self.domain_to_pocket_vector[list(self.domain_to_pocket_vector.keys())[0]])

        x = []
        y = []

        substrate_order = []
        domain_order = []

        counter = 0

        for domain in domains:
            for substrate in self.substrate_to_bitvector:
                substrate_order.append(substrate)
                domain_order.append(domain)

                substrate_order.append(substrate)
                domain_order.append(domain)
                if substrate in self.domain_to_substrates[domain]:
                    interaction = 1
                else:
                    interaction = 0
                datapoint = RFTrainingPoint(domain, substrate, interaction)
                datapoint.build_vector(self.domain_to_pocket_vector, self.substrate_to_bitvector)
                x.append(np.array(datapoint.vector, dtype='uint8'))
                y.append(interaction)
                counter += 1
                if counter % 1000 == 0:
                    print(f"Processed {counter} datapoints.")
        y = np.array(y, dtype='uint8')
        substrate_order = np.array(substrate_order)
        domain_order = np.array(domain_order)

        if under_sample:
            print("Under-sampling..")

            undersampler = RandomUnderSampler(random_state=25051989)

            under_x, under_y = undersampler.fit_resample(x, y)
            under_substrate_order = []
            under_domain_order = []
            for index in undersampler.sample_indices_:
                under_substrate_order.append(substrate_order[index])
                under_domain_order.append(domain_order[index])

            return under_x, under_y, under_substrate_order, under_domain_order

        elif over_sample:
            print("Over-sampling..")
            oversampler = SMOTE(random_state=25051989)
            over_x, over_y = oversampler.fit_resample(x, y)

            over_substrate_order = []
            over_domain_order = []
            for index in oversampler.sample_indices_:
                over_substrate_order.append(substrate_order[index])
                over_domain_order.append(domain_order[index])

            return over_x, over_y, over_substrate_order, over_domain_order

        else:
            return x, y, substrate_order, domain_order


def run():
    args = parse_arguments()

    random_forest = RandomForestParasect(args.b, args.p, args.s, args.d)
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
