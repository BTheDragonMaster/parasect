from sys import argv
import os

from matplotlib import pyplot as plt

from paras.scripts.parsers.parsers import parse_specificities, parse_test_results, parse_substrate_list
from paras.scripts.parsers.iterate_over_dir import iterate_over_dir


def plot_performance_per_aa(specificities, included_substrates, test_results, out_plot):
    plt.figure(figsize=(12, 12))
    substrates = parse_substrate_list(included_substrates)
    domain_to_prediction = parse_test_results(test_results)
    domain_to_spec = parse_specificities(specificities)
    substrate_to_counts = {}
    for substrate in substrates:
        substrate_to_counts[substrate] = {"correct": 0,
                                          "incorrect": 0}

    for domain, prediction in domain_to_prediction.items():
        specs = domain_to_spec[domain]
        if prediction in specs:
            substrate_to_counts[prediction]["correct"] += 1
        else:
            for spec in specs:
                if spec in substrates:
                    substrate_to_counts[spec]["incorrect"] += 1

    substrates = sorted(substrates,
                        key=lambda x: substrate_to_counts[x]["correct"] + substrate_to_counts[x]["incorrect"],
                        reverse=True)

    for i, substrate in enumerate(substrates):
        if i == 0:
            barplot = plt.bar(substrate, [substrate_to_counts[substrate]["correct"],
                                          substrate_to_counts[substrate]["incorrect"]],
                              color=["lightsteelblue", "navajowhite"],
                              bottom=[0, substrate_to_counts[substrate]["correct"]],
                              label=["correct", "incorrect"],
                              width=0.8)
        else:
            barplot = plt.bar(substrate, [substrate_to_counts[substrate]["correct"],
                                          substrate_to_counts[substrate]["incorrect"]],
                              color=["lightsteelblue", "navajowhite"],
                              bottom=[0, substrate_to_counts[substrate]["correct"]],
                              width=0.8)

        labels = [v if v != 0 else "" for v in barplot.datavalues]
        plt.bar_label(barplot, label_type='center', labels=labels)

    plt.xticks(substrates, rotation='vertical')
    plt.legend(loc="upper right", bbox_to_anchor=(1.2, 1.05))
    plt.ylabel("# datapoints")
    plt.savefig(out_plot, bbox_inches='tight')
    plt.clf()


def plot_performance_per_aa_bulk(specificities, included_substrates, results_directory):

    for _, model_dir in iterate_over_dir(results_directory, get_dirs=True):
        for _, extraction_dir in iterate_over_dir(model_dir, get_dirs=True):
            for _, featurisation_dir in iterate_over_dir(extraction_dir, get_dirs=True):
                for test_dir_name, test_dir in iterate_over_dir(featurisation_dir, get_dirs=True):
                    if test_dir_name == 'test_performance':
                        test_results = os.path.join(test_dir, "test_results.txt")
                        out_file = os.path.join(test_dir, "performance_per_substrate.svg")
                        plot_performance_per_aa(specificities, included_substrates, test_results, out_file)




if __name__ == "__main__":
    # plot_performance_per_aa(argv[1], argv[2], argv[3], argv[4])
    plot_performance_per_aa_bulk(argv[1], argv[2], argv[3])