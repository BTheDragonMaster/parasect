from argparse import ArgumentParser
from statistics import mean

from paras.scripts.parsers.parsers import parse_sequence_feature_importances
from paras.scripts.parsers.iterate_over_dir import iterate_over_dir


def parse_arguments():
    parser = ArgumentParser(description="Average multiple feature importance files into one.")
    parser.add_argument('-i', required=True, help="Feature importance file.")
    parser.add_argument('-o', required=True, help="Output file.")
    args = parser.parse_args()

    return args


def run():
    args = parse_arguments()
    all_feature_to_importance = {}
    averaged_feature_to_importance = {}
    for file_name, file_path in iterate_over_dir(args.i, '.txt'):
        feature_to_importance = parse_sequence_feature_importances(file_path)
        for feature, importance in feature_to_importance.items():
            if feature not in all_feature_to_importance:
                all_feature_to_importance[feature] = []
            all_feature_to_importance[feature].append(importance)

    for feature, importances in all_feature_to_importance.items():
        importance = mean(importances)
        averaged_feature_to_importance[feature] = importance

    with open(args.o, 'w') as out:
        out.write('category\timportance\n')
        for feature, importance in averaged_feature_to_importance.items():
            out.write(f"{feature}\t{importance:.5f}\n")


if __name__ == "__main__":
    run()
