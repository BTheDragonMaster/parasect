import os
from typing import Any, Optional
from statistics import mean, stdev

from sklearn.ensemble import RandomForestClassifier
import numpy as np

from sqlalchemy.orm import Session

from parasect.core.chem import fingerprint_to_bitvector
from parasect.database.query_database import get_substrates_from_name
from parasect.database.build_database import AdenylationDomain, Substrate
from parasect.core.parsing import parse_list, parse_pcs
from parasect.core.featurisation import get_domain_features
from parasect.model_training.data_processing.plotting.confusion_matrix import plot_matrix, write_matrix


def write_parasect_results(domain_to_sorted_predictions, included_substrates, out_file):
    with open(out_file, 'w') as out:
        out.write('domain_name')
        for i in range(len(included_substrates)):
            out.write(f"\tprediction_{i}\tinteraction_probability_{i}")
        out.write('\n')
        for domain, sorted_predictions in domain_to_sorted_predictions.items():
            out.write(domain.get_name())
            for prediction in sorted_predictions:
                out.write(f"\t{prediction[1].name}\t{prediction[0]}")
            out.write('\n')


def write_parasect_metrics(tp, tp_probs, fp, fp_probs, tn, tn_probs, fn, fn_probs,
                           correct_ranks, incorrect_ranks, substrate_metrics, out_file):
    with open(out_file, 'w') as out:
        out.write("substrate\tprecision\trecall\tf1\tmean_rank_correct\tsd_rank_correct\tmean_rank_incorrect\tsd_rank_incorrect\ttp\tmean_prob_tp\tsd_prob_tp\tfp\tmean_prob_fp\tsd_prob_fp\tfn\tmean_prob_fn\tsd_prob_fn\ttn\tmean_prob_tn\tsd_prob_tn\n")
        eps = 1e-9



        precision = tp / (tp + fp) if (tp + fp) > eps else "N/A"
        recall = tp / (tp + fn) if (tp + fn) > eps else "N/A"

        print(f"Overall precision: {precision}")
        print(f"Overall recall: {recall}")

        f1 = 2 * (precision * recall) / (precision + recall) if (precision != 'N/A' and recall != 'N/A' and precision + recall > eps) else 'N/A'
        tp_prob, tp_prob_sd = get_mean_stdev(tp_probs)
        fp_prob, fp_prob_sd = get_mean_stdev(fp_probs)
        tn_prob, tn_prob_sd = get_mean_stdev(tn_probs)
        fn_prob, fn_prob_sd = get_mean_stdev(fn_probs)
        rank_correct, rank_correct_sd = get_mean_stdev(correct_ranks)
        rank_incorrect, rank_incorrect_sd = get_mean_stdev(incorrect_ranks)

        out.write(f"all\t{precision}\t{recall}\t{f1}\t{rank_correct}\t{rank_correct_sd}\t{rank_incorrect}\t{rank_incorrect_sd}\t{tp}\t{tp_prob}\t{tp_prob_sd}\t{fp}\t{fp_prob}\t{fp_prob_sd}\t{fn}\t{fn_prob}\t{fn_prob_sd}\t{tn}\t{tn_prob}\t{tn_prob_sd}\n")
        for substrate, metrics in substrate_metrics.items():
            tp = metrics["TP"]
            fp = metrics["FP"]
            tn = metrics["TN"]
            fn = metrics["FN"]

            precision = tp / (tp + fp) if (tp + fp) > eps else "N/A"
            recall = tp / (tp + fn) if (tp + fn) > eps else "N/A"
            f1 = 2 * (precision * recall) / (precision + recall) if (precision != 'N/A' and recall != 'N/A' and precision + recall > eps) else 'N/A'

            tp_prob, tp_prob_sd = get_mean_stdev(metrics["TP probabilities"])
            fp_prob, fp_prob_sd = get_mean_stdev(metrics["FP probabilities"])
            tn_prob, tn_prob_sd = get_mean_stdev(metrics["TN probabilities"])
            fn_prob, fn_prob_sd = get_mean_stdev(metrics["FN probabilities"])

            rank_correct, rank_correct_sd = get_mean_stdev(metrics["correct_ranks"])
            rank_incorrect, rank_incorrect_sd = get_mean_stdev(metrics["incorrect_ranks"])

            out.write(
                f"{substrate.name}\t{precision}\t{recall}\t{f1}\t{rank_correct}\t{rank_correct_sd}\t{rank_incorrect}\t{rank_incorrect_sd}\t{tp}\t{tp_prob}\t{tp_prob_sd}\t{fp}\t{fp_prob}\t{fp_prob_sd}\t{fn}\t{fn_prob}\t{fn_prob_sd}\t{tn}\t{tn_prob}\t{tn_prob_sd}\n")


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


def get_mean_stdev(list_of_values):
    if len(list_of_values) > 0:
        mean_value = mean(list_of_values)
        if len(list_of_values) > 1:
            stdev_value = stdev(list_of_values)
        else:
            stdev_value = "N/A"
    else:
        mean_value = "N/A"
        stdev_value = "N/A"

    return mean_value, stdev_value


def fmt(val):
    return f"{val:.3f}" if isinstance(val, (float, int)) else str(val)


def write_confidences(correct_confidences, incorrect_confidences, per_substrate_confidences, out_file):
    with open(out_file, 'w') as out:
        out.write(f"substrate\tconfidence_correct\tstdev_correct\tconfidence_incorrect\tstdev_incorrect\n")
        mean_correct, stdev_correct = get_mean_stdev(correct_confidences)
        mean_incorrect, stdev_incorrect = get_mean_stdev(incorrect_confidences)

        out.write(
            f"all\t{fmt(mean_correct)}\t{fmt(stdev_correct)}\t{fmt(mean_incorrect)}\t{fmt(stdev_incorrect)}\n"
        )
        substrates = sorted(per_substrate_confidences.keys())
        for substrate in substrates:
            substrate_mean_correct, substrate_stdev_correct = get_mean_stdev(
                per_substrate_confidences[substrate]["correct"])
            substrate_mean_incorrect, substrate_stdev_incorrect = get_mean_stdev(
                per_substrate_confidences[substrate]["incorrect"])

            out.write(
                f"{substrate}\t{fmt(substrate_mean_correct)}\t{fmt(substrate_stdev_correct)}\t"
                f"{fmt(substrate_mean_incorrect)}\t{fmt(substrate_stdev_incorrect)}\n"
            )

def test_parasect_esm(session: Session, classifier: RandomForestClassifier,
                      domains: list[AdenylationDomain],
                      included_substrates_file: str,
                      embeddings_file: str,
                      out_dir: str,
                      hashes: list[int],
                      n_components: int = 100) -> None:
    """Test PARASECT RF trained on ESM PCA

    :param session: database session
    :type session: Session
    :param classifier: PARASECT classifier
    :type classifier: RandomForestClassifier
    :param domains: list of adenylation domains to test
    :type domains: list[AdenylationDomain]
    :param included_substrates_file: path to file containing included substrates
    :type included_substrates_file: str
    :param embeddings_file:
    :type embeddings_file:
    :param out_dir: path to output directory
    :type out_dir: str
    :param hashes: list of hashes representing substrates used for training PARASECT model
    :type hashes: list[int]
    :param n_components:
    :type n_components:
    """
    included_substrate_names = parse_list(included_substrates_file)

    included_substrates = []
    substrates_to_predictions = {}

    domains_pca, pca_data = parse_pcs(embeddings_file)
    domain_to_pcs = {
        domain: pca_data[i][:n_components]
        for i, domain in enumerate(domains_pca)}

    for name in included_substrate_names:
        included_substrates.append(get_substrates_from_name(session, name)[0])

    for substrate in included_substrates:
        features = []
        for domain in domains:
            domain_features = domain_to_pcs[domain.get_name()]
            substrate_features = fingerprint_to_bitvector(hashes, set(substrate.fingerprint))
            features.append(np.array(list(domain_features) + list(substrate_features)))

        features = np.array(features)

        predict_y = classifier.predict_proba(features)

        predictions = []
        for prediction in predict_y:

            prob_per_class = list(zip(prediction, classifier.classes_))
            for prob, cls in prob_per_class:
                if cls == 1:
                    predictions.append(prob)

        substrates_to_predictions[substrate] = predictions

    write_predictions_parasect(domains, substrates_to_predictions, out_dir)


def test_parasect_signatures(session: Session, classifier: RandomForestClassifier,
                             domains: list[AdenylationDomain],
                             included_substrates_file: str,
                             out_dir: str,
                             hashes: list[int]) -> None:
    """Test PARASECT RF trained on signatures

    :param session: database session
    :type session: Session
    :param classifier: PARASECT classifier
    :type classifier: RandomForestClassifier
    :param domains: list of adenylation domains to test
    :type domains: list[AdenylationDomain]
    :param included_substrates_file: path to file containing included substrates
    :type included_substrates_file: str
    :param out_dir: path to output directory
    :type out_dir: str
    :param hashes: list of hashes representing substrates used for training PARASECT model
    :type hashes: list[int]
    """

    included_substrate_names = parse_list(included_substrates_file)

    included_substrates = []
    substrates_to_predictions = {}
    for name in included_substrate_names:
        included_substrates.append(get_substrates_from_name(session, name)[0])

    for substrate in included_substrates:
        features = []
        for domain in domains:

            domain_features = get_domain_features(domain.extended_signature)
            substrate_features = fingerprint_to_bitvector(hashes, set(substrate.fingerprint))
            features.append(np.array(domain_features + substrate_features))

        features = np.array(features)

        predict_y = classifier.predict_proba(features)

        predictions = []
        for prediction in predict_y:

            prob_per_class = list(zip(prediction, classifier.classes_))
            for prob, cls in prob_per_class:
                if cls == 1:
                    predictions.append(prob)

        substrates_to_predictions[substrate] = predictions

    write_predictions_parasect(domains, substrates_to_predictions, out_dir)


def test_paras_signatures(classifier: RandomForestClassifier,
                          domains: list[AdenylationDomain], included_substrates_file: str,
                          out_dir: str) -> None:
    """Test PARAS RF trained on signatures

    :param classifier: random forest classifier
    :type classifier: RandomForestClassifier
    :param domains: list of adenylation domains
    :type domains: list[AdenylationDomain]
    :param included_substrates_file: path to file containing included substrates
    :type included_substrates_file: str
    :param out_dir: path to output directory
    :type out_dir: str
    """

    features = np.array([get_domain_features(domain.extended_signature) for domain in domains])
    predict_y = classifier.predict_proba(features)


    predictions = []
    for prediction in predict_y:
        prob_per_class = list(zip(prediction, classifier.classes_))

        prob_per_class.sort(reverse=True)
        predictions.append(prob_per_class)

    write_predictions_paras(domains, predictions, included_substrates_file, out_dir)


def test_paras_esm(classifier: RandomForestClassifier,
                   domains: list[AdenylationDomain],
                   embeddings_file: str, included_substrates_file: str, out_dir: str,
                   n_components: int = 100) -> None:
    """Test PARAS RF trained on ESM PCA embeddings

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

    domains_pca, pca_data = parse_pcs(embeddings_file)

    name_to_row = {name: row for name, row in zip(domains_pca, pca_data)}
    features = np.array([name_to_row[d.get_name()][:n_components] for d in domains])

    predict_y = classifier.predict_proba(features)

    predictions = []
    for prediction in predict_y:
        prob_per_class = list(zip(prediction, classifier.classes_))

        prob_per_class.sort(reverse=True)
        predictions.append(prob_per_class)

    write_predictions_paras(domains, predictions, included_substrates_file, out_dir)


def write_predictions_paras(domains: list[AdenylationDomain], predictions: list,
                            included_substrates_file: str,
                            out_dir: str):
    """Write predictions to output directory

    :param domains: list of adenylation domains
    :type domains: list[AdenylationDomain]
    :param predictions: predictions made by classifier
    :type predictions: NDArray
    :param included_substrates_file: path to file containing included substrates
    :type included_substrates_file: str
    :param out_dir: path to output directory
    :type out_dir: str
    """
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    print("Writing files...")

    overall_correct = 0
    overall_incorrect = 0
    per_substrate_counts = {}  # {substrate_name: [correct, total]}
    per_substrate_confidences = {}

    included_substrate_names = parse_list(included_substrates_file)

    for substrate in included_substrate_names:
        per_substrate_counts[substrate] = [0, 0]
        per_substrate_confidences[substrate] = {"correct": [],
                                                "incorrect": []}

    substrate_nr = len(included_substrate_names)

    confusion_matrix = np.zeros(shape=(substrate_nr, substrate_nr))

    correct_confidences = []
    incorrect_confidences = []

    for i, prediction_matrix in enumerate(predictions):
        prediction = prediction_matrix[0][1]
        confidence = prediction_matrix[0][0]
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
            correct_confidences.append(confidence)
            per_substrate_confidences[prediction]["correct"].append(confidence)
        else:
            overall_incorrect += 1
            incorrect_confidences.append(confidence)
            per_substrate_confidences[prediction]["incorrect"].append(confidence)

        for substrate in substrate_names:
            if substrate in included_substrate_names:
                per_substrate_counts[substrate][1] += 1  # total
            if prediction == substrate:
                per_substrate_counts[substrate][0] += 1  # correct
                break

    # Compute accuracies
    overall_accuracy = overall_correct / (overall_correct + overall_incorrect)
    print(f"Overall accuracy: {overall_accuracy}")
    per_substrate_accuracy = {
        s: (correct / total if total > 0 else None)
        for s, (correct, total) in per_substrate_counts.items()
    }

    out_plot = os.path.join(out_dir, "confusion_matrix.svg")
    out_matrix = os.path.join(out_dir, "confusion_matrix.txt")
    out_accuracy = os.path.join(out_dir, "accuracy.txt")
    out_confidence = os.path.join(out_dir, "confidence.txt")

    plot_matrix(confusion_matrix, included_substrate_names, out_plot)
    write_matrix(confusion_matrix, included_substrate_names, out_matrix)
    write_confidences(correct_confidences, incorrect_confidences, per_substrate_confidences, out_confidence)
    write_accuracy(overall_accuracy, per_substrate_accuracy, out_accuracy)


def reverse_rank(values: list[float], idx: int) -> int:
    """
    Compute 1-indexed reverse rank of element at index `idx` in `values`.
    Largest value gets rank 1.
    """
    target = values[idx]
    # Count how many elements are strictly greater than target
    rank = sum(v > target for v in values) + 1
    return rank


def write_predictions_parasect(domains: list[AdenylationDomain], substrate_to_predictions: dict[Substrate, list],
                               out_dir: str):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    print("Writing files..")

    substrate_to_metrics: dict[Substrate, dict[str, Any]] = {}
    substrate_nr = len(substrate_to_predictions)
    confusion_matrix = np.zeros(shape=(substrate_nr, substrate_nr))
    domain_to_best_prob: dict[AdenylationDomain, float] = {}
    domain_to_best: dict[AdenylationDomain, Optional[Substrate]] = {}
    domain_to_probs: dict[AdenylationDomain, dict[str, list[float]]] = {}
    domain_to_predictions: dict[AdenylationDomain, list[float]] = {}
    domain_to_sorted_predictions: dict[AdenylationDomain, list[tuple[float, Substrate]]] = {}
    for domain in domains:
        domain_to_best_prob[domain] = 0.0
        domain_to_best[domain] = None
        domain_to_probs[domain] = {"correct": [],
                                   "incorrect": []}
        domain_to_predictions[domain] = []
        domain_to_sorted_predictions[domain] = []

    substrates = sorted(substrate_to_predictions.keys())

    for substrate in substrates:
        predictions = substrate_to_predictions[substrate]
        substrate_to_metrics[substrate] = {"TP": 0,
                                           "FP": 0,
                                           "TN": 0,
                                           "FN": 0,
                                           "TP probabilities": [],
                                           "FP probabilities": [],
                                           "TN probabilities": [],
                                           "FN probabilities": [],
                                           "correct_ranks": [],
                                           "incorrect_ranks": []
                                           }

        for i, interaction_probability in enumerate(predictions):
            domain = domains[i]
            domain_to_predictions[domain].append(interaction_probability)
            domain_to_sorted_predictions[domain].append((interaction_probability, substrate))
            if interaction_probability > domain_to_best_prob[domain]:
                domain_to_best[domain] = substrate
                domain_to_best_prob[domain] = interaction_probability
            if substrate in domain.substrates:
                domain_to_probs[domain]["correct"].append(interaction_probability)
            else:
                domain_to_probs[domain]["incorrect"].append(interaction_probability)

    for domain in domains:
        domain_to_sorted_predictions[domain].sort(reverse=True)

    sorted_out = os.path.join(out_dir, "parasect_predictions.txt")
    write_parasect_results(domain_to_sorted_predictions, substrates, sorted_out)

    correct_ranks = []
    incorrect_ranks = []

    tp = 0
    fp = 0
    fn = 0
    tn = 0

    tp_probs = []
    fp_probs = []
    fn_probs = []
    tn_probs = []

    for i, domain in enumerate(domains):

        for j, substrate in enumerate(substrates):

            rank = reverse_rank(list(domain_to_predictions[domain]), j)
            if substrate in domain.substrates:
                substrate_to_metrics[substrate]["correct_ranks"].append(rank)
                correct_ranks.append(rank)
            else:
                substrate_to_metrics[substrate]["incorrect_ranks"].append(rank)
                incorrect_ranks.append(rank)

        best_prediction = domain_to_best[domain]
        prediction_index = substrates.index(best_prediction)
        substrate_index = None
        for substrate in domain.substrates:
            if substrate in substrates:
                substrate_index = substrates.index(substrate)
                break

        if substrate_index is not None:
            confusion_matrix[prediction_index][substrate_index] += 1
        if best_prediction in domain.substrates:
            substrate_to_metrics[best_prediction]["TP"] += 1
            substrate_to_metrics[best_prediction]["TP probabilities"].append(domain_to_best_prob[domain])
            tp += 1
            tp_probs.append(domain_to_best_prob[domain])
            for substrate in substrate_to_predictions:
                if substrate not in domain.substrates:
                    substrate_to_metrics[substrate]["TN"] += 1
                    interaction_probability = substrate_to_predictions[substrate][i]
                    substrate_to_metrics[best_prediction]["TN probabilities"].append(interaction_probability)
                    tn += 1
                    tn_probs.append(interaction_probability)

        else:
            substrate_to_metrics[best_prediction]["FP"] += 1
            substrate_to_metrics[best_prediction]["FP probabilities"].append(domain_to_best_prob[domain])
            fp += 1
            fp_probs.append(domain_to_best_prob[domain])

            for substrate in domain.substrates:
                if substrate in substrates:
                    substrate_to_metrics[substrate]["FN"] += 1
                    interaction_probability = substrate_to_predictions[substrate][i]
                    substrate_to_metrics[best_prediction]["FN probabilities"].append(interaction_probability)
                    fn += 1
                    fn_probs.append(interaction_probability)

    out_file = os.path.join(out_dir, "parasect_performance.txt")

    write_parasect_metrics(tp, tp_probs, fp, fp_probs, tn, tn_probs, fn, fn_probs, correct_ranks, incorrect_ranks,
                           substrate_to_metrics, out_file)

    out_plot = os.path.join(out_dir, "confusion_matrix.svg")
    out_matrix = os.path.join(out_dir, "confusion_matrix.txt")


    plot_matrix(confusion_matrix, substrates, out_plot)
    write_matrix(confusion_matrix, [s.name for s in substrates], out_matrix)
