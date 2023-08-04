from sys import argv
import os
from matplotlib import pyplot as plt

from paras.scripts.parsers.tabular import Tabular
from paras.scripts.parsers.external_parsers import parse_sandpuma_data
from paras.scripts.parsers.parsers import parse_domain_list, parse_specificities, parse_substrate_list


def plot_histogram(summary_file, out_plot, sequence_only=True):
    summary = Tabular(summary_file, [0])
    predictor_to_metrics = {}
    has_no_call = False
    labels = []
    predictors = []
    for datapoint in summary.data:
        true_substrates = summary.get_value(datapoint, "substrate").split('|')

        for predictor in summary.categories[1:-1]:
            if sequence_only:
                if 'structure' in predictor:
                    continue

            if sequence_only and 'sequence' in predictor:
                label = predictor.split(' ')[0]
            else:
                label = predictor

            if predictor not in predictors:
                labels.append(label)
                predictors.append(predictor)

            if predictor not in predictor_to_metrics:
                predictor_to_metrics[predictor] = {"correct": 0,
                                                   "incorrect": 0,
                                                   "no call": 0}
            predictions = summary.get_value(datapoint, predictor).split('|')
            is_correct = False
            no_call = False
            for prediction in predictions:
                if prediction in true_substrates:
                    is_correct = True
                elif prediction in ['no_call', 'N/A', "no_confident_result", "no_force_needed"]:
                    no_call = True
                    has_no_call = True

            if is_correct:
                predictor_to_metrics[predictor]["correct"] += 1
            elif no_call:
                predictor_to_metrics[predictor]["no call"] += 1
            else:
                predictor_to_metrics[predictor]["incorrect"] += 1

    for i, predictor in enumerate(predictors):

        predictor_name = labels[i]

        correct, incorrect, no_call = predictor_to_metrics[predictor]["correct"], \
                                      predictor_to_metrics[predictor]["incorrect"], \
                                      predictor_to_metrics[predictor]["no call"]
        print(predictor_name, correct, incorrect, no_call, correct / (correct + incorrect), correct / (correct + incorrect + no_call))
        total = correct + incorrect + no_call

        correct = correct / total * 100
        incorrect = incorrect / total * 100
        no_call = no_call / total * 100
        if has_no_call:
            bars = [correct, incorrect, no_call]
            bottom = [0, correct, correct + incorrect]
            label = ["correct", "incorrect", "no call"]
            colors = ["lightsteelblue", "navajowhite", "gainsboro"]
        else:
            if not sequence_only:
                bars = [correct]
                bottom = [0]
                colors = ["lightsteelblue"]
                label = ["correct"]
            else:
                bars = [correct, incorrect]
                bottom = [0, correct]
                label = ["correct", "incorrect"]
                colors = ["lightsteelblue", "navajowhite"]

        if i == 0:
            plt.bar(predictor_name, bars, color=colors, label=label, bottom=bottom)
        else:
            plt.bar(predictor_name, bars, color=colors, bottom=bottom)

    plt.xticks(labels, rotation='vertical')

    if not sequence_only:
        plt.ylabel("Accuracy (%)")
        plt.ylim((0, 100))
    else:
        plt.ylabel("Proportion of datapoints (%)")
        plt.legend(loc="upper right", bbox_to_anchor=(1.2, 1.05))

    plt.savefig(out_plot, bbox_inches='tight')


def benchmark(results_dir, prefix="PARAS"):
    structure_results_phylogeny = os.path.join(results_dir,
                                               "phylogeny_structure/test_performance/test_results.txt")
    sequence_results_phylogeny = os.path.join(results_dir,
                                              "phylogeny_sequence/test_performance/test_results.txt")
    sequence_structure_results_phylogeny = os.path.join(results_dir,
                                                        "phylogeny_sequence_structure/test_performance/test_results.txt")
    structure_results_class = os.path.join(results_dir,
                                           "class_structure/test_performance/test_results.txt")
    sequence_results_class = os.path.join(results_dir,
                                          "class_sequence/test_performance/test_results.txt")
    sequence_structure_results_class = os.path.join(results_dir,
                                                    "class_sequence_structure/test_performance/test_results.txt")

    files = [structure_results_phylogeny,
             sequence_results_phylogeny,
             sequence_structure_results_phylogeny,
             structure_results_class,
             sequence_results_class,
             sequence_structure_results_class]


def benchmark_sandpuma(domain_file, parasect_data, included_substrates, force_call=True):
    domain_to_prediction = parse_sandpuma_data()
    domain_list = parse_domain_list(domain_file)
    domain_to_specificity = parse_specificities(parasect_data)
    substrates = parse_substrate_list(included_substrates)

    nr_correct = 0
    nr_incorrect = 0
    nr_no_call = 0
    substrate_not_included = 0
    missing_substrates = set()

    for domain in domain_list:
        correct = False
        has_included_substrate = False
        no_call = False

        if force_call:
            predictions = domain_to_prediction[domain].ensemble
        else:
            predictions = domain_to_prediction[domain].sandpuma

        for prediction in predictions:
            if prediction in substrates:
                has_included_substrate = True
            if prediction in domain_to_specificity[domain]:
                correct = True
            if prediction == 'no_call':
                no_call = True

        if correct:
            nr_correct += 1
        elif no_call:
            nr_no_call += 1
        elif not has_included_substrate:
            substrate_not_included += 1
            for prediction in predictions:
                missing_substrates.add(prediction)
        else:
            nr_incorrect += 1

    print("=============")
    print("SANDPUMA:")
    print("=============")

    print(f"Correct: {nr_correct}")
    print(f"Incorrect: {nr_incorrect}")
    print(f"No call: {nr_no_call}")
    print(f"Substrate not included in PARAS: {substrate_not_included}")
    print("Missing substrates:")
    print(missing_substrates)

    print('\n')

    accuracy = nr_correct / (nr_correct + nr_incorrect)

    print(f"Accuracy excluding no-call: {accuracy}")

    accuracy = nr_correct / (nr_correct + nr_incorrect + nr_no_call)

    print(f"Accuracy including no-call: {accuracy}")


def benchmark_nrpspredictor(domain_file, parasect_data, included_substrates):
    domain_to_prediction = parse_sandpuma_data()
    domain_list = parse_domain_list(domain_file)
    domain_to_specificity = parse_specificities(parasect_data)
    substrates = parse_substrate_list(included_substrates)

    nr_correct = 0
    nr_incorrect = 0
    nr_no_call = 0
    substrate_not_included = 0
    missing_substrates = set()

    for domain in domain_list:
        correct = False
        has_included_substrate = False
        no_call = False

        predictions = domain_to_prediction[domain].svm

        for prediction in predictions:
            if prediction in substrates:
                has_included_substrate = True
            if prediction in domain_to_specificity[domain]:
                correct = True
            if prediction in ['no_call', 'N/A', "no_confident_result", "no_force_needed"]:
                no_call = True

        if correct:
            nr_correct += 1
        elif no_call:
            nr_no_call += 1
        elif not has_included_substrate:
            substrate_not_included += 1
            for prediction in predictions:
                missing_substrates.add(prediction)
        else:
            nr_incorrect += 1

    print("=============")
    print("NRPSPredictor2:")
    print("=============")
    print(f"Correct: {nr_correct}")
    print(f"Incorrect: {nr_incorrect}")
    print(f"No call: {nr_no_call}")
    print(f"Substrate not included in PARAS: {substrate_not_included}")
    print("Missing substrates:")
    print(missing_substrates)

    print('\n')

    accuracy = nr_correct / (nr_correct + nr_incorrect)

    print(f"Accuracy excluding no-call: {accuracy}")

    accuracy = nr_correct / (nr_correct + nr_incorrect + nr_no_call)

    print(f"Accuracy including no-call: {accuracy}")


if __name__ == "__main__":
    # benchmark_sandpuma(argv[1], argv[2], argv[3], False)
    # benchmark_nrpspredictor(argv[1], argv[2], argv[3])
    plot_histogram(argv[1], argv[2], sequence_only=False)

