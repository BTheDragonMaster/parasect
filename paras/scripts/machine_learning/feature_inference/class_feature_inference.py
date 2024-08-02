from typing import Union
import os
import pylab
from sys import argv

from joblib import load
import dtreeviz
import numpy as np

from paras.scripts.parsers.fasta import read_fasta
from paras.scripts.feature_extraction.sequence_feature_extraction.seq_to_features import PROPERTIES_FILE
from paras.scripts.feature_extraction.sequence_feature_extraction.seq_to_features import get_sequence_features
from paras.scripts.parsers.parsers import parse_amino_acid_properties, parse_morgan_fingerprint, parse_specificities, \
    parse_domain_list
import paras.data.sequence_data.sequences
import paras.data

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['font.size'] = 9.0

SIGNATURES = os.path.join(os.path.dirname(paras.data.sequence_data.sequences.__file__), 'active_site_34_hmm.fasta')
DATASET = os.path.join(os.path.dirname(paras.data.__file__), 'parasect_dataset.txt')
DOMAIN_TO_SPEC = parse_specificities(DATASET)
DOMAIN_LIST = os.path.join(os.path.dirname(paras.data.__file__), 'domain_list_hmm.txt')

SUBSTRATE_TO_GROUP = {"2,3-dihydroxybenzoic acid": "non-amino acid",
                      "2,4-diaminobutyric acid": "polar nitrogen",
                      "2-aminoadipic acid": "polar oxygen",
                      "2-aminoisobutyric acid": "double sidechain",
                      "3,5-dihydroxyphenylglycine": "polar aromatic",
                      "4-hydroxyphenylglycine": "polar aromatic",
                      "D-alanine": "small",
                      "N5-formyl-N5-hydroxyornithine": "polar oxygen",
                      "N5-hydroxyornithine": "polar oxygen",
                      "R-beta-hydroxytyrosine": "polar aromatic",
                      "alanine": "small",
                      "anthranilic acid": "non-amino acid",
                      "arginine": "polar nitrogen",
                      "asparagine": "polar oxygen",
                      "aspartic acid": "polar oxygen",
                      "beta-alanine": "small",
                      "cysteine": "cysteine",
                      "glutamic acid": "polar oxygen",
                      "glutamine": "polar oxygen",
                      "glycine": "small",
                      "histidine": "polar aromatic",
                      "isoleucine": "hydrophobic",
                      "leucine": "hydrophobic",
                      "lysine": "polar nitrogen",
                      "ornithine": "polar nitrogen",
                      "phenylalanine": "aromatic",
                      "pipecolic acid": "polar nitrogen",
                      "proline": "small",
                      "salicylic acid": "non-amino acid",
                      "serine": "polar oxygen",
                      "threonine": "polar oxygen",
                      "tryptophan": "aromatic",
                      "tyrosine": "polar aromatic",
                      "valine": "hydrophobic"
                      }

COLOUR_DICT = {"non-amino acid": "pink",
               "double side chain": "mistyrose",
               "polar oxygen": "maroon",
               "polar nitrogen": "lightskyblue",
               "cysteine": "gold",
               "polar aromatic": "orange",
               "aromatic": "lightgrey",
               "hydrophobic": "grey",
               "small": "darkgrey",
               "double sidechain": "gainsboro"}


def count_residues(position):
    aa_to_count = {}

    id_to_seq = read_fasta(SIGNATURES)
    for _, seq in id_to_seq.items():
        aa = seq[position]
        if aa not in aa_to_count:
            aa_to_count[aa] = 0
        aa_to_count[aa] += 1

    return aa_to_count


def get_substrate_group_counts(counts_per_substrate, classes, substrate_groups):
    substrate_group_counts = [0] * len(substrate_groups)
    for i, value_per_class in enumerate(counts_per_substrate):
        substrate_group = SUBSTRATE_TO_GROUP[classes[i]]
        index = substrate_groups.index(substrate_group)
        substrate_group_counts[index] += value_per_class

    return substrate_group_counts


def plot_pycharts(tree, classifier, out_folder, label):
    out = os.path.join(out_folder, label)
    if not os.path.exists(out):
        os.mkdir(out)
    colour_dict = {}
    cm = pylab.get_cmap('gist_rainbow')
    substrate_groups = list(set(SUBSTRATE_TO_GROUP.values()))
    nr_colours = len(substrate_groups)
    for i in range(nr_colours):
        colour = cm(1. * i / nr_colours)
        substrate = substrate_groups[i]
        colour_dict[substrate] = colour

    categories = get_categories()
    for node_index, feature_index in enumerate(tree.feature):
        if feature_index >= 0:
            values_per_class_left = tree.value[tree.children_left[node_index]][0]
            values_per_class_right = tree.value[tree.children_right[node_index]][0]
            substrate_group_counts_left = get_substrate_group_counts(values_per_class_left, classifier.classes_,
                                                                     substrate_groups)
            substrate_group_counts_right = get_substrate_group_counts(values_per_class_right, classifier.classes_,
                                                                      substrate_groups)

            feature = categories[feature_index]

            file_name_left = os.path.join(out, f"node_{node_index}_less_than_{feature}.svg")
            file_name_right = os.path.join(out, f"node_{node_index}_more_than_{feature}.svg")
            # plot_pychart(list(map(int, values_per_class_left)), classifier.classes_, colour_dict, file_name_left)
            # plot_pychart(list(map(int, values_per_class_right)), classifier.classes_, colour_dict, file_name_right)
            plot_pychart(list(map(int, substrate_group_counts_left)), substrate_groups, COLOUR_DICT, file_name_left)
            plot_pychart(list(map(int, substrate_group_counts_right)), substrate_groups, COLOUR_DICT, file_name_right)


def get_data(classes, domain_list):
    x = []
    y = []
    domains = parse_domain_list(domain_list)
    id_to_seq = read_fasta(SIGNATURES)

    for domain in domains:
        seq = id_to_seq[domain]
        features = get_sequence_features(seq)
        substrate = None
        for substrate_option in DOMAIN_TO_SPEC[domain]:
            if substrate_option in classes:
                substrate = substrate_option
                break
        if not substrate:
            print(DOMAIN_TO_SPEC[domain][0])
        assert substrate
        spec_found = False
        for i, spec in enumerate(classes):
            if spec == substrate:
                y.append(i)
                spec_found = True

        assert spec_found

        x.append(features)

    x = np.array(x)
    y = np.array(y)

    return x, y


def plot_pychart(sizes, labels, colour_dict, out_file):

    plotted_sizes = []
    plotted_labels = []
    colours = []

    for i, size in enumerate(sizes):
        if size >= 1:
            plotted_sizes.append(size)
            colours.append(colour_dict[labels[i]])
            plotted_labels.append(labels[i])

    plt.pie(plotted_sizes, colors=colours, labels=plotted_labels)
    plt.tight_layout(pad=0.0)
    plt.savefig(out_file, bbox_inches=0, transparent=True)
    plt.clf()


def make_piechart_tree(classifier, estimator, out_dir, label):
    estimator.classes_ = list(map(int, estimator.classes_))
    categories = get_categories()
    x, y = get_data(classifier.classes_, DOMAIN_LIST)
    viz_model = dtreeviz.model(estimator, x, y, feature_names=categories, target_name="substrate", class_names=classifier.classes_)
    view = viz_model.view(orientation="LR")
    view.save(os.path.join(out_dir, f"{label}_pychart_tree.svg"))


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


CATEGORY_TO_VERBOSE = {"WOLS870101": "Hydrophilicity z1",
                       "WOLS870102": "Size z2",
                       "WOLS870103": "Electron state z3",
                       "FAUJ880109": "Hydrogen donors",
                       "GRAR740102": "Polarity 1",
                       "RADA880108": "Hydrophobicity",
                       "ZIMJ680103": "Polarity 2",
                       "TSAJ990101": "Volume",
                       "CHOP780201": "a-helix",
                       "CHOP780202": "b-sheet",
                       "CHOP780203": "b-turn",
                       "ZIMJ680104": "Isoelectric point",
                       "NEU1": "Hydrophilicity pc1",
                       "NEU2": "Hydrophobicity pc2",
                       "NEU3": "Hydrophobicity pc3"}


def get_categories():
    _, sequence_categories = parse_amino_acid_properties(PROPERTIES_FILE, return_categories=True)
    categories = []

    for i in range(34):
        for category in sequence_categories:
            categories.append(f"{CATEGORY_TO_VERBOSE[category]}|res{i + 1}")
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


def analyse_first_splits(model, out_folder):
    categories = get_categories()
    classifier = load(model)
    feature_to_count = {}
    res_to_count = {}
    seen_features = set()
    for i, estimator in enumerate(classifier.estimators_):
        tree = estimator.tree_
        first_feature = categories[tree.feature[0]]

        if not os.path.exists(out_folder):
            os.mkdir(out_folder)

        res = first_feature.split('|')[-1]
        if res not in res_to_count:
            res_to_count[res] = 0
        if first_feature not in feature_to_count:
            feature_to_count[first_feature] = 0
        feature_to_count[first_feature] += 1
        res_to_count[res] += 1

        if first_feature not in seen_features and feature_to_count[first_feature] > 20:
            seen_features.add(first_feature)
            plot_pycharts(tree, classifier, out_folder, label=first_feature)
            # make_piechart_tree(classifier, estimator, out_folder, label=first_feature)

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
        # print("threonine")
        # pprint(count_residues_per_substrate(position, "threonine"))
        # print("other")
        # pprint(count_residues_per_substrate(position, "threonine", True))


if __name__ == "__main__":
    # get_splits_per_class_paras(PARAS)
    analyse_first_splits(argv[1], argv[2])
