import os
from pprint import pprint
from sys import argv

from statistics import mean

from paras.scripts.parsers.iterate_over_dir import iterate_over_dir
from paras.scripts.parsers.parsers import parse_metrics


def analyse_performance(crossval_folder):
    performance = {}
    for output_type, type_path in iterate_over_dir(crossval_folder, get_dirs=True):
        performance[output_type] = {}
        for split_type, split_path in iterate_over_dir(type_path, get_dirs=True):
            performance[output_type][split_type] = {}
            for feature_type, feature_path in iterate_over_dir(split_path, get_dirs=True):
                performance[output_type][split_type][feature_type] = {"train": {},
                                                                      "test": {}}

                for dir_name, crossval_path in iterate_over_dir(feature_path, get_dirs=True):
                    if 'crossval' in dir_name:
                        test_metrics_file = os.path.join(crossval_path, "test_performance/accuracy.txt")
                        train_metrics_file = os.path.join(crossval_path, "train_performance/accuracy.txt")
                        test_metrics = parse_metrics(test_metrics_file)
                        train_metrics = parse_metrics(train_metrics_file)

                        for metric, value in test_metrics.items():
                            if metric not in performance[output_type][split_type][feature_type]["test"]:
                                performance[output_type][split_type][feature_type]["test"][metric] = []

                            if value != 'N/A':
                                performance[output_type][split_type][feature_type]["test"][metric].append(value)

                        for metric, value in train_metrics.items():
                            if metric not in performance[output_type][split_type][feature_type]["train"]:
                                performance[output_type][split_type][feature_type]["train"][metric] = []
                            if value != 'N/A':
                                performance[output_type][split_type][feature_type]["train"][metric].append(value)


    for output_type, rest in performance.items():
        for split_type, rest1 in rest.items():
            for feature_type, rest2 in rest1.items():
                for set_type, metric_to_values in rest2.items():
                    for metric, values in metric_to_values.items():
                        if values:
                            performance[output_type][split_type][feature_type][set_type][metric] = mean(values)

    pprint(performance)


if __name__ == "__main__":
    analyse_performance(argv[1])
