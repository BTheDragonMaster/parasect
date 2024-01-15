from typing import Union
import os
from pprint import pprint

from joblib import load

from paras.scripts.general import PARAS, PARASECT
from paras.scripts.parsers.fasta import read_fasta
from paras.scripts.feature_extraction.sequence_feature_extraction.seq_to_features import PROPERTIES_FILE
from paras.scripts.parsers.parsers import parse_amino_acid_properties, parse_morgan_fingerprint, parse_specificities
import paras.data.sequence_data.sequences
import paras.data

SIGNATURES = os.path.join(os.path.dirname(paras.data.sequence_data.sequences.__file__), 'active_site_34_hmm.fasta')
DATASET = os.path.join(os.path.dirname(paras.data.__file__), 'trustworthy_parasect_dataset.txt')
DOMAIN_TO_SPEC = parse_specificities(DATASET)


def count_residues(position):
    aa_to_count = {}

    id_to_seq = read_fasta(SIGNATURES)
    for _, seq in id_to_seq.items():
        aa = seq[position]
        if aa not in aa_to_count:
            aa_to_count[aa] = 0
        aa_to_count[aa] += 1

    return aa_to_count


def count_residues_per_substrate(position, substrate, return_counts_other=False):
    aa_to_count_substrate = {}
    aa_to_count_other = {}

    id_to_seq = read_fasta(SIGNATURES)
    for domain, substrates in DOMAIN_TO_SPEC.items():
        if domain in id_to_seq:
            seq = id_to_seq[domain]
            aa = seq[position]

            if substrate in substrates:
                dictionary = aa_to_count_substrate
            else:
                dictionary = aa_to_count_other

            if aa not in dictionary:
                dictionary[aa] = 0
            dictionary[aa] += 1

    if return_counts_other:

        return aa_to_count_other
    else:
        return aa_to_count_substrate


def get_categories():
    _, sequence_categories = parse_amino_acid_properties(PROPERTIES_FILE, return_categories=True)
    categories = []

    for i in range(34):
        for category in sequence_categories:
            categories.append(f"{category}_res{i + 1}")
    return categories


def get_splits_per_class_paras(model):
    categories = get_categories()
    classifier = load(model)

    feature_frequencies = [0.0] * len(categories)

    node_count_matrix: list[list[Union[float, int]]] = []

    for i in range(len(classifier.classes_)):
        row = [0.0] * len(categories)
        node_count_matrix.append(row)

    for tree in [estimator.tree_ for estimator in classifier.estimators_]:
        for node_index, feature_index in enumerate(tree.feature):
            if feature_index >= 0:
                feature_frequencies[feature_index] += 1
                values_per_class = tree.value[node_index][0]
                for j, value in enumerate(values_per_class):
                    if value > 0:
                        node_count_matrix[j][feature_index] += 1

    for i, substrate in enumerate(classifier.classes_):
        relative_importances = node_count_matrix[i]
        normalised_importances = []
        for j, frequency in enumerate(relative_importances):
            if frequency > 0.0:
                # print(categories[j], frequency, feature_frequencies[j])
                normalised_importances.append(frequency / feature_frequencies[j])
            else:
                normalised_importances.append(0.0)

        print(f"{substrate} (normalized)")
        named_features = list(zip(normalised_importances, categories))
        named_raw_features = list(zip(relative_importances, categories))
        named_raw_features.sort(key=lambda x: x[0], reverse=True)
        named_features.sort(key=lambda x: x[0], reverse=True)
        for importance, named_feature in named_features[:10]:
            print(f"{named_feature}: {importance}")

        print(f"{substrate} (raw)")
        for importance, named_feature in named_raw_features[:10]:
            print(f"{named_feature}: {importance}")


def analyse_first_splits(model):
    categories = get_categories()
    classifier = load(model)
    feature_to_count = {}
    res_to_count = {}
    for tree in [estimator.tree_ for estimator in classifier.estimators_]:
        first_feature = categories[tree.feature[0]]
        values_per_class_left = tree.value[tree.children_left[0]][0]
        values_per_class_right = tree.value[tree.children_right[0]][0]

        if first_feature == "CHOP780201_res25" or first_feature == "WOLS870102_res25":
            print('\n')
            print(f"{first_feature} < {tree.threshold[0]}")
            for i, substrate in enumerate(classifier.classes_):
                print(f"{substrate}\t{values_per_class_left[i]}\t{values_per_class_right[i]}")


        res = first_feature.split('_')[-1]
        if res not in res_to_count:
            res_to_count[res] = 0
        if first_feature not in feature_to_count:
            feature_to_count[first_feature] = 0
        feature_to_count[first_feature] += 1
        res_to_count[res] += 1

    features_counts = []
    for feature, count in feature_to_count.items():
        features_counts.append((feature, count))

    res_counts = []
    for res, count in res_to_count.items():
        res_counts.append((res, count))

    features_counts.sort(key=lambda x: x[1], reverse=True)
    res_counts.sort(key=lambda x: x[1], reverse=True)

    for feature, count in features_counts:
        print(f"{feature}\t{count}")

    for res, count in res_counts:
        print(f"{res}\t{count}")
        position = int(res.split('res')[-1]) - 1
        print("threonine")
        pprint(count_residues_per_substrate(position, "threonine"))
        print("other")
        pprint(count_residues_per_substrate(position, "threonine", True))




def get_splits_per_class_parasect(model, fingerprint_file):
    _, sequence_categories = parse_amino_acid_properties(PROPERTIES_FILE, return_categories=True)
    _, compound_categories = parse_morgan_fingerprint(fingerprint_file, return_categories=True)


if __name__ == "__main__":
    # get_splits_per_class_paras(PARAS)
    analyse_first_splits(PARAS)