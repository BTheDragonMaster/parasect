from sys import argv
import os

import matplotlib.pyplot as plt
import numpy as np

from paras.scripts.parsers.tabular import Tabular
from paras.scripts.data_analysis.plotting.plot_confusion_matrices import plot_matrix

SIGNATURE_POS = [210, 213, 214, 230, 234, 235, 236, 237, 238, 239, 240, 243, 278, 279, 299, 300, 301, 302,
                 303, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334]
STACH_POS = [235, 236, 239, 278, 299, 301, 322, 330, 331]


def compile_extraction_data(extraction_file_1, extraction_file_2, out_file):
    domains = []
    extraction_data_1 = Tabular(extraction_file_1, [0])

    for datapoint in extraction_data_1.data:
        domain = extraction_data_1.get_value(datapoint, "Domain")
        domains.append(domain)

    with open(out_file, 'w') as out:
        with open(extraction_file_1, 'r') as extraction_1:
            for line in extraction_1:
                out.write(line)
        with open(extraction_file_2, 'r') as extraction_2:
            extraction_2.readline()
            for line in extraction_2:
                domain = line.split('\t')[0]
                if domain not in domains:
                    out.write(line)


def compare_extraction_methods(extraction_file, out_dir, mode='stach'):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    extraction_data = Tabular(extraction_file, [0])

    nr_mismatches = 0
    nr_total_mismatches = 0
    nr_domains = 0

    residues = "-ACDEFGHIKLMNPQRSTVWY"
    cm = []
    position_to_mismatch = {}

    if mode == 'stach':
        for position in STACH_POS:
            position_to_mismatch[position] = 0
    elif mode == 'signature':
        for position in SIGNATURE_POS:
            position_to_mismatch[position] = 0
    else:
        raise Exception("Mode has to be stach or signature.")

    sequence_gap_counts = 0
    structure_gap_counts = 0

    for _ in residues:
        row = []
        for _ in residues:
            row.append(0)

        cm.append(row)

    for datapoint in extraction_data.data:
        nr_domains += 1

        sequence_gaps = int(extraction_data.get_value(datapoint, "gap_count_sequence_based"))
        sequence_gap_counts += sequence_gaps
        structure_gaps = int(extraction_data.get_value(datapoint, "gap_count_structure_based"))
        structure_gap_counts += structure_gaps

        mismatching_positions = extraction_data.get_value(datapoint, "mismatching_positions").split(', ')
        for position in mismatching_positions:
            if position:
                position_to_mismatch[int(position) + 66] += 1

        mismatches = extraction_data.get_value(datapoint, "mismatches").split(', ')
        for mismatch in mismatches:
            if mismatch:
                nr_total_mismatches += 1

                sequence_res = mismatch[0]
                structure_res = mismatch[1]
                if sequence_res != '-' and structure_res != '-':
                    nr_mismatches += 1
                cm[residues.index(sequence_res)][residues.index(structure_res)] += 1

    out_matrix = os.path.join(out_dir, "mismatch_matrix.svg")

    plot_matrix(cm, list(residues), out_matrix, cluster=False)

    positions = sorted(list(position_to_mismatch.keys()))
    counts = []

    for position in positions:
        counts.append(position_to_mismatch[position])

    out_bar = os.path.join(out_dir, "mismatches_per_position.svg")

    plt.bar(np.arange(len(positions)), counts, color='lightsteelblue')
    plt.xticks(np.arange(len(positions)), positions, rotation='vertical')
    plt.xlabel("Active site residue")
    plt.ylabel("# mismatches")
    plt.savefig(out_bar, bbox_inches='tight')

    average_mismatches = nr_mismatches / nr_domains
    average_mismatches_gaps = nr_total_mismatches / nr_domains
    average_gaps_sequence = sequence_gap_counts / nr_domains
    average_gaps_structure = structure_gap_counts / nr_domains
    print("Average number of mismatches:", average_mismatches)
    print("Average number of mismatches including gaps:", average_mismatches_gaps)
    print("Average number of gaps in sequence:", average_gaps_sequence)
    print("Average number of gaps in structure:", average_gaps_structure)


if __name__ == "__main__":
    compare_extraction_methods(argv[1], argv[2], argv[3])
