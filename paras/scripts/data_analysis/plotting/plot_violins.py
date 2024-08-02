import os
from sys import argv

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import argparse

import paras.data.train_test_splits
import paras.data.benchmarking
import paras.data.benchmarking.result_summaries

from paras.scripts.parsers.parsers import parse_domain_list, parse_closest_identities
from paras.scripts.parsers.tabular import Tabular


benchmarking_dir = os.path.dirname(paras.data.benchmarking.__file__)

SUMMARY_FILE_CLASS = os.path.join(os.path.dirname(paras.data.benchmarking.result_summaries.__file__), "class.txt")
SUMMARY_FILE_PHYLOGENY = os.path.join(os.path.dirname(paras.data.benchmarking.result_summaries.__file__), "phylogeny.txt")

DOMAINS_CLASS_SANDPUMA = os.path.join(benchmarking_dir, 'benchmark_class_filtered.txt')
DOMAINS_PHYLOGENY_SANDPUMA = os.path.join(benchmarking_dir, 'benchmark_phylogeny_filtered.txt')

CLASS_IDENTITIES = os.path.join(os.path.dirname(paras.data.train_test_splits.__file__), 'class/domain_to_identity.txt')
PHYLOGENY_IDENTITIES = os.path.join(os.path.dirname(paras.data.train_test_splits.__file__), 'phylogeny/domain_to_identity.txt')
SANDPUMA_IDENTITIES = os.path.join(benchmarking_dir, 'benchmarking_set/pid.res.tsv')
NRPSPREDICTOR_PHYLOGENY_IDENTITIES = os.path.join(benchmarking_dir, 'benchmarking_set/pid_nrpspredictor_phylogeny.txt')
NRPSPREDICTOR_CLASS_IDENTITIES = os.path.join(benchmarking_dir, 'benchmarking_set/pid_nrpspredictor_class.txt')


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', type=str, help="Path to PARAS & PARASECT results folder.")
    parser.add_argument('-o', type=str, help="Path to output file.")
    parser.add_argument('-d', type=str, help="Path to domain list representing holdout test set.")
    parser.add_argument('-m', type=str, default='class', help="Method used to create holdout test set. Must be 'class' or 'phylogeny'.")
    parser.add_argument('-s', action='store_true', help="Add sandpuma results to summary file.")
    parser.add_argument('-a', action='store_true', help="Add adenpredictor results to summary file.")
    parser.add_argument('-n', action='store_true', help="Add nrpspredictor results to summary file.")

    args = parser.parse_args()

    return args


def plot_violins(identity_file, nrpspredictor_identity_file, summary_file, domain_list, out_file):
    included_domains = parse_domain_list(domain_list)
    summary = Tabular(summary_file, [0])
    domain_to_identity = parse_closest_identities(identity_file)
    domain_to_sandpuma_identity = parse_closest_identities(SANDPUMA_IDENTITIES)
    domain_to_nrpspredictor_identity = parse_closest_identities(nrpspredictor_identity_file)

    metrics = {"prediction": [],
               "closest % identity": [],
               "domain": [],
               "predictor": []}

    for datapoint in summary.data:
        domain = summary.get_value(datapoint, "domain")
        if domain in included_domains:

            true_substrates = summary.get_value(datapoint, "substrate").split('|')

            paras_predictors = []
            other_predictors = []

            for predictor in summary.categories[2:]:
                if 'Structure' in predictor:
                    continue
                elif 'SANDPUMA' in predictor and (predictor != "SANDPUMA (ensemble)"): #and predictor != "SANDPUMA"):
                    continue
                elif 'Sequence' in predictor:
                    paras_predictors.append(predictor)
                else:
                    other_predictors.append(predictor)

            paras_predictors.sort(key=lambda x: len(x))
            other_predictors.sort()

            predictors = paras_predictors + other_predictors

            for predictor in predictors:

                label = predictor

                if "PARAS" in predictor or "PARASECT" in predictor:
                    identity = domain_to_identity[domain]
                    label = predictor.split(' ')[0]
                elif predictor == "NRPSPredictor2":
                    identity = domain_to_nrpspredictor_identity[domain]
                else:
                    if "SANDPUMA" in predictor:
                        label = 'SANDPUMA'
                    identity = domain_to_sandpuma_identity[domain]

                metrics["domain"].append(domain)
                metrics["predictor"].append(label)
                metrics["closest % identity"].append(identity)

                predictions = summary.get_value(datapoint, predictor).split('|')
                is_correct = False
                no_call = False

                for prediction in predictions:
                    if prediction in true_substrates:
                        is_correct = True
                    elif prediction in ['no_call', 'N/A', "no_confident_result", "no_force_needed"]:
                        no_call = True

                if is_correct:
                    metrics["prediction"].append("correct")
                elif no_call:
                    # pass
                    metrics["prediction"].append("incorrect")

                else:
                    metrics["prediction"].append("incorrect")

    plt.figure(figsize=(4, 3.2))
    palette = {"correct": "steelblue",
               "incorrect": "lightcoral",
               "no call": "gainsboro"}

    df = pd.DataFrame(data=metrics)
    ax = sns.violinplot(data=df, x="predictor", y="closest % identity", hue="prediction", scale="count", cut=0,
                        linewidth=1.0, palette=palette, split=True, bw=0.2)
    plt.xticks(rotation=90)
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    ax.set(xlabel=None)
    plt.savefig(out_file, bbox_inches='tight')

    plt.clf()


if __name__ == "__main__":
    out_dir = argv[1]
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    class_out = os.path.join(out_dir, "class_violin.svg")
    phylogeny_out = os.path.join(out_dir, "phylogeny_violin.svg")
    plot_violins(CLASS_IDENTITIES, NRPSPREDICTOR_CLASS_IDENTITIES, argv[2], argv[3], class_out)
    plot_violins(PHYLOGENY_IDENTITIES, NRPSPREDICTOR_PHYLOGENY_IDENTITIES, argv[4], argv[5], phylogeny_out)
