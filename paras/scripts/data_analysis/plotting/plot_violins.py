import os
from sys import argv

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

import paras.data.train_test_splits
import paras.data.benchmarking
import paras.data.benchmarking.result_summaries

from paras.scripts.parsers.parsers import parse_domain_list, parse_closest_identities
from paras.scripts.parsers.tabular import Tabular


benchmarking_dir = os.path.dirname(paras.data.benchmarking.__file__)

DOMAINS_CLASS = os.path.join(os.path.dirname(paras.data.train_test_splits.__file__), 'class/test_filtered.txt')
DOMAINS_PHYLOGENY = os.path.join(os.path.dirname(paras.data.train_test_splits.__file__), 'phylogeny/test_filtered.txt')
SUMMARY_FILE_CLASS = os.path.join(os.path.dirname(paras.data.benchmarking.result_summaries.__file__), "class.txt")
SUMMARY_FILE_PHYLOGENY = os.path.join(os.path.dirname(paras.data.benchmarking.result_summaries.__file__), "phylogeny.txt")

DOMAINS_CLASS_SANDPUMA = os.path.join(benchmarking_dir, 'benchmark_class_filtered.txt')
DOMAINS_PHYLOGENY_SANDPUMA = os.path.join(benchmarking_dir, 'benchmark_phylogeny_filtered.txt')

CLASS_IDENTITIES = os.path.join(os.path.dirname(paras.data.train_test_splits.__file__), 'class/domain_to_identity.txt')
PHYLOGENY_IDENTITIES = os.path.join(os.path.dirname(paras.data.train_test_splits.__file__), 'phylogeny/domain_to_identity.txt')
SANDPUMA_IDENTITIES = os.path.join(benchmarking_dir, 'benchmarking_set/pid.res.tsv')


def plot_violins(identity_file, summary_file, domain_list, out_file):
    included_domains = parse_domain_list(domain_list)
    summary = Tabular(summary_file, [0])
    domain_to_identity = parse_closest_identities(identity_file)
    domain_to_sandpuma_identity = parse_closest_identities(SANDPUMA_IDENTITIES)

    metrics = {"prediction": [],
               "closest % identity": [],
               "domain": [],
               "predictor": []}

    for datapoint in summary.data:
        domain = summary.get_value(datapoint, "domain")
        if domain in included_domains:

            true_substrates = summary.get_value(datapoint, "substrate").split('|')

            for predictor in summary.categories[1:-1]:
                if 'structure' in predictor:
                    continue
                elif 'SANDPUMA' in predictor and (predictor != "SANDPUMA (ensemble)"): #and predictor != "SANDPUMA"):
                    continue
                elif 'sequence' in predictor:
                    label = predictor.split(' ')[0]
                else:
                    label = predictor

                if "PARAS" in predictor or "PARASECT" in predictor:
                    identity = domain_to_identity[domain]
                else:
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
                    pass
                    # metrics["prediction"].append("incorrect")

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
    plot_violins(CLASS_IDENTITIES, SUMMARY_FILE_CLASS, DOMAINS_CLASS_SANDPUMA, class_out)
    plot_violins(PHYLOGENY_IDENTITIES, SUMMARY_FILE_PHYLOGENY, DOMAINS_PHYLOGENY_SANDPUMA, phylogeny_out)
