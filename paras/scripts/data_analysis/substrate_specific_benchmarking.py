import os
from sys import argv
import paras.data.compound_data
from paras.scripts.parsers.parsers import parse_summary_file, parse_substrate_list
from paras.scripts.parsers.iterate_over_dir import iterate_over_dir
from matplotlib import pyplot as plt
from paras.scripts.parsers.tabular import Tabular
import numpy as np

SUBSTRATE_FILE = os.path.join(os.path.dirname(paras.data.compound_data.__file__), 'included_substrates.txt')


def write_substrate_metrics(summary_file, included_substrates, output_folder):
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    predictor_to_performance = parse_summary_file(summary_file)

    for predictor, performance in predictor_to_performance.items():
        performance.set_per_substrate_metrics(included_substrates)
        output_path = os.path.join(output_folder, f"{predictor}.txt")
        with open(output_path, 'w') as out:
            out.write(f"substrate\ttp\tfp\tfn\tprecision\trecall\tf1-score\n")
            for substrate in included_substrates:
                metrics = performance.substrate_metrics[substrate]
                out.write(f"{substrate}\t{metrics.tp}\t{metrics.fp}\t{metrics.fn}\t{metrics.precision}\t{metrics.recall}\t{metrics.f1}\n")

def make_plot(labels, title0, title1, column0, column1, out_file):
    index = []
    for i, label in enumerate(labels):
        index.append(i + 1)

    hfont = {'fontname': 'Verdana'}
    color0 = 'lightsteelblue'
    color1 = 'navajowhite'

    fig, axes = plt.subplots(figsize=(10, 15), ncols=2, sharey=True)
    fig.tight_layout()

    axes[0].barh(index, column0, align='center', color=color0)
    axes[0].set_title(title0, fontsize=18, pad=15, color='black', **hfont)

    axes[1].barh(index, column1, align='center', color=color1)
    axes[1].set_title(title1, fontsize=18, pad=15, color='black', **hfont)

    axes[0].invert_xaxis()
    plt.gca().invert_yaxis()

    axes[0].set(yticks=index, yticklabels=labels)
    axes[0].yaxis.tick_left()

    axes[0].tick_params(axis='y', colors='black')

    axes[0].xaxis.set_major_locator(plt.MaxNLocator(5))
    axes[1].xaxis.set_major_locator(plt.MaxNLocator(5))
    axes[0].xaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))
    axes[1].xaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))

    # axes[0].set_xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    #
    # axes[1].set_xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

    plt.subplots_adjust(wspace=0, top=0.85, bottom=0.1, left=0.25, right=0.95)

    plt.savefig(out_file)
    plt.clf()

def plot_horizontal(labels, title0, title1, column0, column1, out_file):
    index = []
    for i, label in enumerate(labels):
        index.append(i + 1)

    hfont = {'fontname': 'Arial'}
    color0 = 'lightsteelblue'
    color1 = 'navajowhite'

    fig, ax = plt.subplots(layout='constrained')
    x = np.arange(5)

    substrates = ["L-Thr", "L-Ser", "L-Lys", "L-Phe", "L-Trp"]
    new_column0 = [column0[labels.index("threonine")],
                   column0[labels.index("serine")],
                   column0[labels.index("lysine")],
                   column0[labels.index("phenylalanine")],
                   column0[labels.index("tryptophan")]]
    new_column1 = [column1[labels.index("threonine")],
                   column1[labels.index("serine")],
                   column1[labels.index("lysine")],
                   column1[labels.index("phenylalanine")],
                   column1[labels.index("tryptophan")]]

    columns = [new_column0, new_column1]
    labels = [title0, title1]
    colors = [color0, color1]

    for i, column in enumerate(columns):
        print(column)

        offset = i * 0.25
        rects = ax.bar(x + offset, column, 0.25, label=labels[i], color=colors[i])
        ax.bar_label(rects, padding=3)

    ax.set_xticks(x + 0.25, substrates)
    ax.legend(loc='upper left', ncols=2)
    ax.set_ylim(0, 1)

    plt.savefig(out_file)


def plot_substrate_metrics_horizontal(summary_folder, output_folder, main_set="PARAS (Sequence)"):
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    main_data = None
    for predictor_name, data_path in iterate_over_dir(summary_folder, '.txt'):
        if predictor_name == main_set:
            main_data = Tabular(data_path, [0])

    assert main_data

    for predictor_name, data_path in iterate_over_dir(summary_folder, '.txt'):
        if predictor_name != main_set:
            recall_plots = os.path.join(output_folder, 'recall')
            f1_plots = os.path.join(output_folder, 'f1')
            precision_plots = os.path.join(output_folder, 'precision')
            for folder in [recall_plots, f1_plots, precision_plots]:
                if not os.path.exists(folder):
                    os.mkdir(folder)

            data = Tabular(data_path, [0])

            if 'Structure' not in predictor_name:
                title0 = 'PARAS'
                if 'Sequence' in predictor_name:
                    title1 = predictor_name.split(" (Sequence")[0]
                else:
                    title1 = predictor_name
            else:
                title0 = main_set
                title1 = predictor_name

            out_precision = os.path.join(precision_plots, f"{title0}_{title1}.svg")
            out_recall = os.path.join(recall_plots, f"{title0}_{title1}.svg")
            out_f1 = os.path.join(f1_plots, f"{title0}_{title1}.svg")

            precision0 = list(map(float, main_data.get_column('precision')))
            precision1 = list(map(float, data.get_column('precision')))

            recall0 = list(map(float, main_data.get_column('recall')))
            recall1 = list(map(float, data.get_column('recall')))

            f10 = list(map(float, main_data.get_column('f1-score')))
            f11 = list(map(float, data.get_column('f1-score')))

            labels = main_data.get_column('substrate')
            plot_horizontal(labels, title0, title1, f10, f11, out_f1)


def plot_substrate_metrics(summary_folder, output_folder, main_set="PARAS (Sequence)"):
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    main_data = None
    for predictor_name, data_path in iterate_over_dir(summary_folder, '.txt'):
        if predictor_name == main_set:
            main_data = Tabular(data_path, [0])

    assert main_data

    for predictor_name, data_path in iterate_over_dir(summary_folder, '.txt'):
        if predictor_name != main_set:
            recall_plots = os.path.join(output_folder, 'recall')
            f1_plots = os.path.join(output_folder, 'f1')
            precision_plots = os.path.join(output_folder, 'precision')
            for folder in [recall_plots, f1_plots, precision_plots]:
                if not os.path.exists(folder):
                    os.mkdir(folder)

            data = Tabular(data_path, [0])

            if 'Structure' not in predictor_name:
                title0 = 'PARAS'
                if 'Sequence' in predictor_name:
                    title1 = predictor_name.split(" (Sequence")[0]
                else:
                    title1 = predictor_name
            else:
                title0 = main_set
                title1 = predictor_name

            out_precision = os.path.join(precision_plots, f"{title0}_{title1}.svg")
            out_recall = os.path.join(recall_plots, f"{title0}_{title1}.svg")
            out_f1 = os.path.join(f1_plots, f"{title0}_{title1}.svg")

            precision0 = list(map(float, main_data.get_column('precision')))
            precision1 = list(map(float, data.get_column('precision')))

            recall0 = list(map(float, main_data.get_column('recall')))
            recall1 = list(map(float, data.get_column('recall')))

            f10 = list(map(float, main_data.get_column('f1-score')))
            f11 = list(map(float, data.get_column('f1-score')))

            labels = main_data.get_column('substrate')

            make_plot(labels, title0, title1, precision0, precision1, out_precision)
            make_plot(labels, title0, title1, f10, f11, out_f1)
            make_plot(labels, title0, title1, recall0, recall1, out_recall)



if __name__ == "__main__":
    # write_substrate_metrics(argv[1], parse_substrate_list(SUBSTRATE_FILE), argv[2])
    plot_substrate_metrics_horizontal(argv[1], argv[2])