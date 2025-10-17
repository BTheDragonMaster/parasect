from sys import argv
import os
from argparse import ArgumentParser, Namespace

from sqlalchemy import create_engine
from sqlalchemy.orm import Session

from parasect.core.tabular import Tabular
from parasect.database.query_database import get_domains_from_synonym
import matplotlib.pyplot as plt
import math

def parse_arguments() -> Namespace:
    parser = ArgumentParser()
    parser.add_argument('-db', type=str, required=True, help="Path to PARASECT db")
    parser.add_argument('-bac', type=str, required=True, help="Path to bacterial PARASECT results")
    parser.add_argument('-all', type=str, required=True, help="Path to all PARASECT results")
    parser.add_argument('-out', type=str, required=True, help="Path to output file")
    parser.add_argument('-n', type=int, default=5, help="Top range")
    args = parser.parse_args()
    return args


def count_correct_in_top_x(session, test_results: str, top_x: int):
    test_data = Tabular(test_results)
    domain_to_interactions = {}

    correct = 0
    incorrect = 0

    for datapoint in test_data.rows:
        domain_name = test_data.get_row_value(datapoint, "domain_name")
        domain = get_domains_from_synonym(session, domain_name.split('|')[0])[0]
        top_predictions = [test_data.get_row_value(datapoint, f"prediction_{i}") for i in range(top_x)]

        substrates = [s.name for s in domain.substrates]
        if set(top_predictions).intersection(substrates):
            correct += 1
        else:
            incorrect += 1


    accuracy = correct / (correct + incorrect)

    print(f"Frequency of substrates occurring in {top_x}:")

    print(f"Correct: {correct}")
    print(f"Incorrect: {incorrect}")
    print(f"Accuracy: {accuracy}")

    return accuracy


def plot_performance(session, test_file, top_range):
    x = []
    y = []

    for i in range(top_range):
        accuracy = count_correct_in_top_x(session, test_file, i + 1)
        x.append(i + 1)
        y.append(accuracy * 100)

    plt.plot(x, y)
    plt.scatter(x, y)
    for index in range(len(x)):
        plt.text(x[index] + 0.1, y[index] - 1.5, f"{y[index]:.2f}", size=12,
                 fontdict=dict(size=10),
                 bbox=dict(facecolor="lightsteelblue", alpha=0.5))
    plt.xlabel("Nr predictions considered")
    plt.ylabel("% true substrate in predictions")
    plt.xticks([1, 2, 3, 4, 5], ['1', '2', '3', '4', '5'])
    plt.xlim(1, max(x) + 1)

    nearest_ten = math.floor(min(y) / 10.0) * 10
    if min(y) - nearest_ten < 3:
        nearest_ten -= 5
    plt.ylim(nearest_ten, 100)
    yticks = []

    if nearest_ten % 10 == 0:
        ytick = nearest_ten
    else:
        ytick = nearest_ten + 5

    while ytick <= 100:
        yticks.append(ytick)
        ytick += 10

    yticks_labels = list(map(str, yticks))

    plt.yticks(yticks, yticks_labels)


def plot_performance_bulk(session, all_results, bacterial_results, out_file, top_range=5):
    plot_performance(session, all_results, top_range)
    plot_performance(session, bacterial_results, top_range)
    plt.savefig(out_file, bbox_inches='tight')
    plt.clf()


if __name__ == "__main__":
    args = parse_arguments()
    engine = create_engine(f"sqlite:///{args.db}")
    with Session(engine) as session:
        plot_performance_bulk(session, args.all, args.bac, args.out)

