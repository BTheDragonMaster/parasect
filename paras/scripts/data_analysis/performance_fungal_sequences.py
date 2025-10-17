import os
from sys import argv
from matplotlib import pyplot as plt
from paras.scripts.parsers.parsers import parse_substrate_list, parse_taxonomy
from paras.scripts.parsers.tabular import Tabular

import paras.data.compound_data
import paras.data

SUBSTRATE_FILE = os.path.join(os.path.dirname(paras.data.compound_data.__file__), 'all_substrates.txt')
TAXONOMY_FILE = os.path.join(os.path.dirname(paras.data.__file__), 'taxonomy.txt')


def plot_histogram_fungal(summary_file, out_plot, sequence_only=True, paras_only=False):
    summary = Tabular(summary_file, [0])
    predictor_to_metrics = {}
    has_no_call = False
    labels = []
    predictors = []
    substrates = parse_substrate_list(SUBSTRATE_FILE)
    protein_to_taxonomy = parse_taxonomy(TAXONOMY_FILE)

    for datapoint in summary.data:
        domain_name = summary.get_value(datapoint, "domain")
        protein_name = '.'.join(domain_name.split('.')[:-1])
        if protein_name in protein_to_taxonomy and 'Fungi' in protein_to_taxonomy[protein_name]:
            true_substrates = summary.get_value(datapoint, "substrate").split('|')

            for predictor in summary.categories[2:]:
                if sequence_only:
                    if 'Structure' in predictor:
                        continue

                if paras_only:
                    if 'PARAS' not in predictor:
                        continue

                if sequence_only and 'Sequence' in predictor:
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
                    elif prediction not in substrates:
                        print(predictor, prediction)

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
        # print(predictor_name, correct, incorrect, no_call, correct / (correct + incorrect), correct / (correct + incorrect + no_call))
        total = correct + incorrect + no_call
        print(total)

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

if __name__ == "__main__":
    plot_histogram_fungal(argv[1], argv[2], sequence_only=True, paras_only=False)