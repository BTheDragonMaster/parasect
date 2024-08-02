import os
from sys import argv
import paras.data.compound_data
from paras.scripts.parsers.parsers import parse_summary_file, parse_substrate_list

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


def plot_substrate_metrics(summary_file, included_substrates, output_folder):
    pass


if __name__ == "__main__":
    write_substrate_metrics(argv[1], parse_substrate_list(SUBSTRATE_FILE), argv[2])