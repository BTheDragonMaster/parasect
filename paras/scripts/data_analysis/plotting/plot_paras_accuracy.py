import os
from sys import argv

import matplotlib.pyplot as plt
import numpy as np
from numpy import argmax, sqrt


from paras.scripts.parsers.parsers import parse_specificities, parse_test_results
from paras.scripts.parsers.iterate_over_dir import iterate_over_dir


def plot_per_threshold(specificities, out_dir, step_size, test_results):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    domain_to_prediction = parse_test_results(test_results, get_confidence=True)
    # domain_to_prediction_2 = parse_test_results(test_results_2, get_confidence=True)
    domain_to_spec = parse_specificities(specificities)

    threshold = 0.0

    tprs = []
    fprs = []
    fnrs = []
    precisions = []
    recalls = []
    accuracies = []
    thresholds = []
    nocalls = []

    while threshold <= 1.0:
        thresholds.append(threshold)
        TP = 0
        FP = 0
        FN = 0
        TN = 0

        correct = 0
        incorrect = 0
        nocall = 0
        for domain, prediction in domain_to_prediction.items():
            predicted_spec, confidence = prediction
            true_specs = domain_to_spec[domain]
            if confidence >= threshold:
                if predicted_spec in true_specs:
                    TP += 1
                    correct += 1
                else:
                    incorrect += 1
                    FP += 1
            else:
                nocall += 1
                if predicted_spec in true_specs:
                    FN += 1
                else:
                    TN += 1

        # for domain, prediction in domain_to_prediction_2.items():
        #     predicted_spec, confidence = prediction
        #     true_specs = domain_to_spec[domain]
        #     if confidence >= threshold:
        #         if predicted_spec in true_specs:
        #             TP += 1
        #             correct += 1
        #         else:
        #             incorrect += 1
        #             FP += 1
        #     else:
        #         nocall += 1
        #         if predicted_spec in true_specs:
        #             FN += 1
        #         else:
        #             TN += 1

        tpr = TP / (TP + FN)
        fpr = FP / (FP + TN)

        fnrs.append(FN / (FN + TP))

        nocalls.append(nocall / (correct + incorrect + nocall))

        accuracy = correct / (correct + incorrect)
        accuracies.append(accuracy)
        tprs.append(tpr)
        fprs.append(fpr)
        precisions.append(TP / (TP + FP))
        recalls.append(TP / (TP + FN))
        threshold += step_size

    recalls = np.array(recalls)
    precisions = np.array(precisions)
    tprs = np.array(tprs)
    fprs = np.array(fprs)

    fscore = (2 * precisions * recalls) / (precisions + recalls)
    index = argmax(fscore)

    plt.plot(recalls, precisions)
    plt.scatter(recalls[index], precisions[index], marker='o', color='black', label='Recommended threshold')
    plt.text(recalls[index] - 0.65, precisions[index] - 0.05, f"Recommended threshold: {thresholds[index]:.2f}",
             size=12,
             fontdict=dict(size=10),
             bbox=dict(facecolor="lightsteelblue", alpha=0.5)
             )
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.ylim([0.0, 1.0])
    plt.xlim([0.0, 1.0])
    plt.savefig(os.path.join(out_dir, "pr_curve.svg"), bbox_inches='tight')
    plt.clf()

    gmeans = sqrt(tprs * (1 - fprs))
    ix = argmax(gmeans)
    plt.plot(fprs, tprs)
    plt.scatter(fprs[ix], tprs[ix], marker='o', color='black', label='Recommended threshold')
    plt.text(fprs[ix] + 0.05, tprs[ix] - 0.05, f"Recommended threshold: {thresholds[ix]:.2f}",
             size=12,
             fontdict=dict(size=10),
             bbox=dict(facecolor="lightsteelblue", alpha=0.5)
             )
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.ylim([0.0, 1.0])
    plt.xlim([0.0, 1.0])
    plt.savefig(os.path.join(out_dir, "roc_curve.svg"), bbox_inches='tight')
    plt.clf()

    accuracies = np.array(accuracies)
    thresholds = np.array(thresholds)
    nocalls = np.array(nocalls)

    plt.plot(thresholds, accuracies * 100, label="Accuracy")
    plt.plot(thresholds, nocalls * 100, label="No call")
    # plt.plot(thresholds, fnrs, label="False negative rate")
    # plt.plot(thresholds, tprs, label="True positive rate")

    plt.xlabel('Threshold')
    plt.ylabel('% of data')
    plt.ylim([0, 100])
    plt.xlim([0.0, max(thresholds)])
    plt.legend()
    plt.savefig(os.path.join(out_dir, "accuracy_curve.svg"), bbox_inches='tight')
    plt.clf()


def bin_by_confidence(specificities, out_dir, n_bins, test_results):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    bin_size = 1.0 / n_bins

    domain_to_prediction = parse_test_results(test_results, get_confidence=True)
    # domain_to_prediction_2 = parse_test_results(test_results_2, get_confidence=True)
    domain_to_spec = parse_specificities(specificities)

    bin_ranges = []
    lower_bound = 0.0
    bins = []

    for i in range(n_bins):
        upper_bound = lower_bound + bin_size
        bin_ranges.append((round(lower_bound, 1), round(upper_bound, 1)))
        lower_bound = upper_bound
        bins.append([])

    print(bin_ranges)

    for domain, prediction in domain_to_prediction.items():
        predicted_spec, confidence = prediction
        true_specs = domain_to_spec[domain]
        correct = False
        if predicted_spec in true_specs:
            correct = True
        for i, bin_range in enumerate(bin_ranges):
            if bin_range[0] < confidence <= bin_range[1]:
                bins[i].append(correct)
                break
    #
    # for domain, prediction in domain_to_prediction_2.items():
    #     predicted_spec, confidence = prediction
    #     true_specs = domain_to_spec[domain]
    #     correct = False
    #     if predicted_spec in true_specs:
    #         correct = True
    #     for i, bin_range in enumerate(bin_ranges):
    #         if bin_range[0] < confidence <= bin_range[1]:
    #             bins[i].append(correct)
    #             break

    out_bar = os.path.join(out_dir, 'bin_bar_plot.svg')
    accuracies = []

    for i, data_bin in enumerate(bins[:]):
        if len(data_bin) == 0:
            del bins[i]
            del bin_ranges[i]
        else:
            accuracy = data_bin.count(True) / len(data_bin)
            accuracies.append(accuracy)

    accuracies = np.array(accuracies)

    bin_names = [f"{x[0]}-{x[1]}" for x in bin_ranges]
    positions = np.arange(len(bin_ranges))

    plt.bar(positions, accuracies * 100, color='lightsteelblue')

    plt.xticks(positions, bin_names, rotation='vertical')
    plt.xlabel("Confidence")
    plt.ylabel("Accuracy (%)")
    plt.savefig(out_bar, bbox_inches='tight')
    plt.clf()


def plot_paras_accuracy_bulk(results_directory, specificities, bins=10, step_size=0.01):
    for model_name, model_dir in iterate_over_dir(results_directory, get_dirs=True):
        if model_name == 'paras':
            for _, extraction_dir in iterate_over_dir(model_dir, get_dirs=True):
                for _, featurisation_dir in iterate_over_dir(extraction_dir, get_dirs=True):
                    for test_dir_name, test_dir in iterate_over_dir(featurisation_dir, get_dirs=True):
                        if test_dir_name == 'test_performance':
                            out_dir = os.path.join(test_dir, 'plots')
                            if not os.path.exists(out_dir):
                                os.mkdir(out_dir)
                            test_results = os.path.join(test_dir, 'test_results.txt')
                            plot_per_threshold(specificities, out_dir, step_size, test_results)
                            bin_by_confidence(specificities, out_dir, bins, test_results)



if __name__ == "__main__":
    # plot_per_threshold(argv[1], argv[2], float(argv[3]), argv[4], argv[6])
    # bin_by_confidence(argv[1], argv[2], int(argv[5]), argv[4], argv[6])
    plot_paras_accuracy_bulk(argv[1], argv[2])
