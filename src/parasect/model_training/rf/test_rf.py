import os

from sklearn.ensemble import RandomForestClassifier
import numpy as np
from numpy.typing import NDArray

from parasect.database.build_database import AdenylationDomain
from parasect.core.parsing import parse_list, parse_pcs
from parasect.core.featurisation import get_domain_features
from parasect.model_training.data_processing.plotting.confusion_matrix import plot_matrix, write_matrix


def write_accuracy(overall_accuracy, per_substrate_accuracy, out_file):
    with open(out_file, 'w') as out:
        out.write(f"substrate\taccuracy\n")
        out.write(f"all\t{overall_accuracy:.3f}\n")
        substrates = sorted(per_substrate_accuracy.keys())
        for substrate in substrates:
            accuracy = per_substrate_accuracy[substrate]
            if accuracy is not None:
                out.write(f"{substrate}\t{accuracy:.3f}\n")
            else:
                out.write(f"{substrate}\tN/A\n")


def test_paras_signatures(classifier: RandomForestClassifier,
                          domains: list[AdenylationDomain], included_substrates_file: str,
                          out_dir: str) -> None:
    """Test RF trained on PARAS signatures

    :param classifier: random forest classifier
    :type classifier: RandomForestClassifier
    :param domains: list of adenylation domains
    :type domains: list[AdenylationDomain]
    :param included_substrates_file: path to file containing included substrates
    :type included_substrates_file: str
    :param out_dir: path to output directory
    :type out_dir: str
    """
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    features = np.array([get_domain_features(domain.extended_signature) for domain in domains])
    predictions = classifier.predict(features)

    write_predictions(domains, predictions, included_substrates_file, out_dir)


def test_paras_esm(classifier: RandomForestClassifier,
                   domains: list[AdenylationDomain],
                   embeddings_file: str, included_substrates_file: str, out_dir: str) -> None:
    """Test RF trained on PARAS signatures

    :param classifier: random forest classifier
    :type classifier: RandomForestClassifier
    :param domains: list of adenylation domains
    :type domains: list[AdenylationDomain]
    :param embeddings_file: path to file containing esm embeddings
    :type embeddings_file: str
    :param included_substrates_file: path to file containing included substrates
    :type included_substrates_file: str
    :param out_dir: path to output directory
    :type out_dir: str
    """
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    domains_pca, pca_data = parse_pcs(embeddings_file)

    name_to_row = {name: row for name, row in zip(domains_pca, pca_data)}
    features = np.array([name_to_row[d.get_name()] for d in domains])

    predictions = classifier.predict(features)
    write_predictions(domains, predictions, included_substrates_file, out_dir)


def write_predictions(domains: list[AdenylationDomain], predictions: NDArray[np.str_],
                      included_substrates_file: str,
                      out_dir: str):
    """Write predictions to output directory

    :param domains: list of adenylation domains
    :type domains: list[AdenylationDomain]
    :param predictions: predictions made by classifier
    :type predictions: NDArray[np.str_]
    :param included_substrates_file: path to file containing included substrates
    :type included_substrates_file: str
    :param out_dir: path to output directory
    :type out_dir: str
    """
    overall_correct = 0
    overall_incorrect = 0
    per_substrate_counts = {}  # {substrate_name: [correct, total]}

    included_substrate_names = parse_list(included_substrates_file)

    for substrate in included_substrate_names:
        per_substrate_counts[substrate] = [0, 0]

    substrate_nr = len(included_substrate_names)

    confusion_matrix = np.zeros(shape=(substrate_nr, substrate_nr))

    for i, prediction in enumerate(predictions):
        prediction_index = included_substrate_names.index(prediction)
        substrate_names = [s.name for s in domains[i].substrates]
        substrate_index = None
        for substrate in substrate_names:
            if substrate in included_substrate_names:
                substrate_index = included_substrate_names.index(substrate)
                break

        if substrate_index is not None:
            confusion_matrix[prediction_index][substrate_index] += 1

        if prediction in substrate_names:
            overall_correct += 1
        else:
            overall_incorrect += 1

        for substrate in substrate_names:
            if substrate in included_substrate_names:
                per_substrate_counts[substrate][1] += 1  # total
            if prediction == substrate:
                per_substrate_counts[substrate][0] += 1  # correct

    # Compute accuracies
    overall_accuracy = overall_correct / (overall_correct + overall_incorrect)
    print(f"Overall accuracy: {overall_accuracy}")
    per_substrate_accuracy = {
        s: (correct / total if total > 0 else None)
        for s, (correct, total) in per_substrate_counts.items()
    }

    out_plot = os.path.join(out_dir, "confusion_matrix.svg")
    out_matrix = os.path.join(out_dir, "confusion_matrix.txt")
    out_accuracy = os.path.join(out_dir, "metrics.txt")

    plot_matrix(confusion_matrix, included_substrate_names, out_plot)
    write_matrix(confusion_matrix, included_substrate_names, out_matrix)

    write_accuracy(overall_accuracy, per_substrate_accuracy, out_accuracy)
