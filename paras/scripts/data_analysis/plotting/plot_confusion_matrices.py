import os
from sys import argv
from collections import OrderedDict

import scipy.cluster.hierarchy as hc
import pandas as pd
import scipy.spatial as sp
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from pprint import pprint

from paras.scripts.parsers.parsers import parse_cm_matrix


def write_matrix(matrix, labels, out_file):
    with open(out_file, 'w') as out:
        out.write('Predicted\\Actual\t')
        out.write('\t'.join(labels))
        out.write("\n")
        for i, aa in enumerate(labels):
            out.write(f"{aa}")
            for j, aa_2 in enumerate(labels):
                out.write(f"\t{matrix[i][j]}")
            out.write("\n")


def plot_matrix(matrix, labels, out_file, cluster=True):
    data = pd.DataFrame(matrix, index=labels, columns=labels)

    data_correlation = data.T.corr()
    # distance_matrix = 1 - data_correlation

    # linkage = hc.linkage(sp.distance.squareform(distance_matrix, checks=False), method='average')

    if cluster:

        clustermap = sns.clustermap(data, cmap=sns.cm.rocket_r) #row_linkage=linkage, col_linkage=linkage)
        clustermap.savefig(out_file)
    else:
        clustermap = sns.heatmap(data, cmap=sns.cm.rocket_r)
        plt.savefig(out_file)

    plt.clf()
    plt.close()


def get_clustering(in_file):
    matrix, labels = parse_cm_matrix(in_file)
    data = pd.DataFrame(matrix, index=labels, columns=labels)
    print(data.size)

    data_correlation = data.T.corr()
    print(data_correlation)
    distance_matrix = 1 - data_correlation
    distance_matrix = distance_matrix.fillna(0)

    linkage = hc.linkage(sp.distance.squareform(distance_matrix), method='average')
    return linkage


def plot_cm_matrix(in_file, linkage, out_file):
    matrix, labels = parse_cm_matrix(in_file)
    data = pd.DataFrame(matrix, index=labels, columns=labels)

    clustermap = sns.clustermap(data, cmap=sns.cm.rocket_r, row_linkage=linkage, col_linkage=linkage)

    clustermap.savefig(out_file)
    plt.clf()
    plt.close()


def get_frame_paths(input_folder):
    frame_paths = []
    for file_name in os.listdir(input_folder):
        file_path = os.path.join(input_folder, file_name)
        if file_name.endswith("cm.txt") and os.path.isfile(file_path):
            epoch = int(file_name.split('_')[1])
            frame_paths.append((epoch, file_path))

    frame_paths.sort(key=lambda x: x[0])

    return frame_paths


def make_cm_frames(input_folder, output_folder):
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    frame_paths = get_frame_paths(input_folder)

    last_frame = frame_paths[-1][1]
    linkage = get_clustering(last_frame)

    for epoch, frame in frame_paths:
        frame_out = os.path.join(output_folder, f"{epoch:04d}.png")
        plot_cm_matrix(frame, linkage, frame_out)


def get_averaged_matrix_per_epoch(frames):
    epoch_to_matrices = OrderedDict()
    for frame_set in frames:
        for epoch, frame in frame_set:
            if epoch not in epoch_to_matrices:
                epoch_to_matrices[epoch] = []
            matrix, labels = parse_cm_matrix(frame)
            epoch_to_matrices[epoch].append(matrix)

    epoch_to_matrix = OrderedDict()

    for epoch, matrices in epoch_to_matrices.items():
        average_matrix = average_matrices(matrices)
        epoch_to_matrix[epoch] = average_matrix

    return epoch_to_matrix


def average_matrices(matrices):
    nr_matrices = len(matrices)
    averaged_matrix = matrices[0]
    nr_rows = 0
    nr_columns = 0
    for matrix in matrices[1:]:
        nr_rows = len(matrix)
        for i, row in enumerate(matrix):
            nr_columns = len(row)
            for j, value in enumerate(row):
                averaged_matrix[i][j] += value

    for i in range(nr_rows):
        for j in range(nr_columns):
            averaged_matrix[i][j] = averaged_matrix[i][j] / nr_matrices

    return averaged_matrix


def plot_all_models(models_dir):
    for model_name in os.listdir(models_dir):
        print(model_name)
        model_dir = os.path.join(models_dir, model_name)
        if os.path.isdir(model_dir):
            for extraction_method in os.listdir(model_dir):
                print(extraction_method)
                extraction_dir = os.path.join(model_dir, extraction_method)
                if os.path.isdir(extraction_dir):
                    for featurisation_method in os.listdir(extraction_dir):
                        print(featurisation_method)

                        crossval_dir = os.path.join(extraction_dir, featurisation_method)
                        if os.path.isdir(crossval_dir):
                            out_path = os.path.join(crossval_dir, "confusion_matrix.svg")
                            matrices = []
                            test_matrix, labels = parse_cm_matrix(os.path.join(crossval_dir, "test_performance/confusion_matrix.txt"))

                            for crossval_set in os.listdir(crossval_dir):
                                print(crossval_set)
                                if 'crossval' in crossval_set:
                                    crossval_set_dir = os.path.join(crossval_dir, crossval_set)

                                    if os.path.isdir(crossval_set_dir):

                                        confusion_matrix = os.path.join(crossval_set_dir, "test_performance/confusion_matrix.txt")

                                        matrix, _ = parse_cm_matrix(confusion_matrix)
                                        matrices.append(matrix)

                            average_matrix = average_matrices(matrices)

                            data = pd.DataFrame(average_matrix, index=labels, columns=labels)

                            data_correlation = data.T.corr()

                            data_correlation = data_correlation.fillna(0)
                            distance_matrix = 1 - data_correlation
                            np.fill_diagonal(distance_matrix.values, 0.0)

                            linkage = hc.linkage(sp.distance.squareform(distance_matrix), method='average')

                            test_data = pd.DataFrame(test_matrix, index=labels, columns=labels)

                            clustermap = sns.clustermap(test_data, cmap=sns.cm.rocket_r, row_linkage=linkage,
                                                        col_linkage=linkage, figsize=(10, 10),
                                                        cbar_pos=(0, .6, .03, .4),
                                                        dendrogram_ratio=0.15)

                            clustermap.ax_row_dendrogram.set_visible(False)

                            clustermap.savefig(out_path)
                            plt.clf()
                            plt.close()


def plot_average(crossval_dir, output_folder):
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    last_frames = []
    frames = []
    for directory in os.listdir(crossval_dir):
        directory_path = os.path.join(crossval_dir, directory)
        if os.path.isdir(directory_path) and 'crossval' in directory:
            crossval_set = int(directory.split('_')[1])
            aa_metrics_path = os.path.join(directory_path, "amino_acid_metrics")
            if os.path.isdir(aa_metrics_path):
                frame_paths = get_frame_paths(aa_metrics_path)
                last_frame = frame_paths[-1][1]
                last_frames.append(last_frame)
                frames.append(frame_paths)

    epoch_to_matrix = get_averaged_matrix_per_epoch(frames)
    last_frame_matrices = []
    labels = []

    for frame in last_frames:
        matrix, labels = parse_cm_matrix(frame)
        last_frame_matrices.append(matrix)

    last_frame_matrix = average_matrices(last_frame_matrices)

    data = pd.DataFrame(last_frame_matrix, index=labels, columns=labels)

    data_correlation = data.T.corr()
    distance_matrix = 1 - data_correlation

    linkage = hc.linkage(sp.distance.squareform(distance_matrix), method='average')

    for epoch, matrix in epoch_to_matrix.items():
        frame_out = os.path.join(output_folder, f"{epoch:04d}.png")
        data = pd.DataFrame(matrix, index=labels, columns=labels)

        clustermap = sns.clustermap(data, cmap=sns.cm.rocket_r, row_linkage=linkage, col_linkage=linkage)

        clustermap.savefig(frame_out)
        plt.clf()
        plt.close()


if __name__ == "__main__":
    # make_cm_frames(argv[1], argv[2])
    # plot_average(argv[1], argv[2])
    # linkage = get_clustering(argv[1])
    # plot_cm_matrix(argv[1], linkage, argv[2])
    plot_all_models(argv[1])
