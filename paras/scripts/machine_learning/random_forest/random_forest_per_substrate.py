from paras.scripts.parsers.parsers import parse_moment_vectors, parse_interactions, parse_domain_list
from paras.scripts.data_processing.get_domain_list import make_domain_mapping
from joblib import dump, load
import os

from sklearn.ensemble import RandomForestClassifier
from imblearn.under_sampling import RandomUnderSampler
from imblearn.over_sampling import SMOTE
import argparse
from statistics import mean, median
from typing import Optional


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, required=True, help="Path to input directory containing training and test data.")
    parser.add_argument('-m', type=str, required=True, help="Path to moment file.")
    parser.add_argument('-d', type=str, required=True, help="Path to domain file")
    parser.add_argument('-a', type=str, required=True, help="Path to directory containing interaction files per substrate")
    parser.add_argument('-o', type=str, required=True, help="Path to output directory")
    parser.add_argument('-n', type=int, default=1000, help="Nr of trees in RF")
    parser.add_argument('-t', type=int, default=None, help="Nr of threads. -1 uses all available cores")
    parser.add_argument('-sampling', type=str, default=None, help="balanced, balanced_subsample, over_sample, under_sample or None")

    args = parser.parse_args()
    return args


class RandomForestParasectSingle:
    def __init__(self, substrate, interaction_file, crossval_nr):
        print(interaction_file)
        self.substrate = substrate
        self.crossval_nr = crossval_nr
        self.domain_to_interaction = parse_interactions(interaction_file)

    def train_model(self, train_domain_file, out_path, domain_to_moment_vector, sampling_method, n_jobs: Optional[int] = None,
                    n_estimators: int = 1000, oob_score: bool = True):
        print("Building data..")

        if sampling_method == 'over_sample':
            train_x, train_y = self.build_data(train_domain_file, domain_to_moment_vector, over_sample=True)
        elif sampling_method == 'under_sample':
            train_x, train_y = self.build_data(train_domain_file, domain_to_moment_vector, under_sample=True)
        else:
            train_x, train_y = self.build_data(train_domain_file, domain_to_moment_vector)

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

    def test_model(self, test_domain_file, classifier_file, domain_to_moment_vector):
        test_x, test_y = self.build_data(test_domain_file, domain_to_moment_vector)
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

        print(f"Total points for substrate {self.substrate}: {len(test_y)}")

        if probs_tp:
            print(f"True positive: {TP} (average probability of interaction: {mean(probs_tp)}; median: {median(probs_tp)})")
        else:
            print(f"True positive: 0")
        if probs_fp:
            print(f"False positive: {FP} (average probability of interaction: {mean(probs_fp)}; median: {median(probs_fp)})")
        else:
            print(f"False positive: 0")
        if probs_tn:
            print(f"True negative: {TN} (average probability of interaction: {mean(probs_tn)}; median: {median(probs_tn)})")
        else:
            print(f"True negative: 0")
        if probs_fn:
            print(f"False negative: {FN} (average probability of interaction: {mean(probs_fn)}; median: {median(probs_fn)})")
        else:
            print(f"False negative: 0")
        print('\n')

    def build_data(self, domain_file: str, domain_to_moment_vector, under_sample=False, over_sample=False):
        domains = parse_domain_list(domain_file)
        if under_sample:
            assert not over_sample
        if over_sample:
            assert not under_sample

        x = []
        y = []

        for domain in domains:
            interaction = self.domain_to_interaction[domain]
            moment_vector = domain_to_moment_vector[domain]
            x.append(moment_vector)
            y.append(interaction)

        if under_sample:
            print("Under-sampling..")

            undersampler = RandomUnderSampler(random_state=25051989)

            under_x, under_y = undersampler.fit_resample(x, y)

            return under_x, under_y

        elif over_sample:
            print("Over-sampling..")
            oversampler = SMOTE(random_state=25051989)
            over_x, over_y = oversampler.fit_resample(x, y)

            return over_x, over_y

        else:
            return x, y


def run():
    args = parse_arguments()

    domain_mapping = make_domain_mapping(args.d)
    old_domain_to_moment_vector = parse_moment_vectors(args.m)
    domain_to_moment_vector = {}
    for old_enzyme, vector in old_domain_to_moment_vector.items():
        enzyme = domain_mapping[old_enzyme]
        domain_to_moment_vector[enzyme] = vector

    if not os.path.exists(args.o):
        os.mkdir(args.o)

    for substrate_name in os.listdir(args.i):
        substrate_dir = os.path.join(args.i, substrate_name)
        if os.path.isdir(substrate_dir):
            crossvalidation_dir = os.path.join(substrate_dir, "crossvalidation")
            if os.path.exists(crossvalidation_dir) and os.path.isdir(crossvalidation_dir):
                interaction_file = os.path.join(args.a, f"{substrate_name}.txt")
                assert os.path.exists(interaction_file) and os.path.isfile(interaction_file)
                nr_to_train_test_pairs = {}
                for domain_file in os.listdir(crossvalidation_dir):
                    if domain_file.endswith('.txt'):
                        domain_file_name = domain_file.split('.txt')[0]
                        _, data_type, crossval_nr = domain_file_name.split('_')
                        crossval_nr = int(crossval_nr)
                        if crossval_nr not in nr_to_train_test_pairs:
                            nr_to_train_test_pairs[crossval_nr] = {}
                            nr_to_train_test_pairs[crossval_nr]["train"] = None
                            nr_to_train_test_pairs[crossval_nr]["test"] = None
                        nr_to_train_test_pairs[crossval_nr][data_type] = os.path.join(crossvalidation_dir, domain_file)

                for crossval_set in nr_to_train_test_pairs:
                    train_domains = nr_to_train_test_pairs[crossval_set]["train"]
                    test_domains = nr_to_train_test_pairs[crossval_set]["test"]
                    assert train_domains is not None and test_domains is not None

                    classifier_out = os.path.join(args.o, f"{substrate_name}_{crossval_set}.parasect")

                    random_forest = RandomForestParasectSingle(substrate_name, interaction_file, crossval_set)
                    random_forest.train_model(train_domains, classifier_out, domain_to_moment_vector, sampling_method=args.sampling,
                                              n_estimators=args.n, n_jobs=args.t)
                    random_forest.test_model(test_domains, classifier_out, domain_to_moment_vector)


if __name__ == "__main__":
    run()
