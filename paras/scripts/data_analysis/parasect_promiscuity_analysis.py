from dataclasses import dataclass
from sys import argv

from paras.scripts.parsers.parsers import parse_specificities, parse_domain_list, parse_substrate_list
from paras.scripts.parsers.tabular import Tabular

import matplotlib
import matplotlib.pyplot as plt
import numpy as np


def autolabel(rects, ax):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')

def plot_promiscuity(out_file):

    labels = ['substrate-stratified', 'taxonomy-stratified']
    correct = [13, 10]
    partially_correct = [4, 11]
    incorrect = [1, 1]

    x = np.arange(len(labels))  # the label locations
    width = 0.2  # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(x - width, correct, width, label='Correct', color='lightsteelblue')
    rects2 = ax.bar(x, partially_correct, width, label='Partially correct', color='gainsboro')
    rects3 = ax.bar(x + width, incorrect, width, label='Incorrect', color='navajowhite')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('#Datapoints')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45)
    ax.legend(loc='upper right', bbox_to_anchor=(1.0, 1.3))

    autolabel(rects1, ax)
    autolabel(rects2, ax)
    autolabel(rects3, ax)
    ax.set_ylim(0, 13)

    fig.tight_layout()

    plt.savefig(out_file)
    plt.clf()


def get_promiscuous(domain_list, parasect_data, included_substrate_list=None):
    domains = parse_domain_list(domain_list)
    domain_to_spec = parse_specificities(parasect_data)
    filtered_domain_to_spec = {}
    included_substrates = parse_substrate_list(included_substrate_list)

    for domain in list(domain_to_spec.keys()):
        spec = domain_to_spec[domain]
        new_spec = []
        for substrate in spec:
            if substrate in included_substrates:
                new_spec.append(substrate)
        if new_spec:
            filtered_domain_to_spec[domain] = new_spec

    promiscuous_domain_to_spec = {}
    for domain in domains:
        if domain in filtered_domain_to_spec and len(filtered_domain_to_spec[domain]) > 1:
            promiscuous_domain_to_spec[domain] = filtered_domain_to_spec[domain]
    return promiscuous_domain_to_spec


@dataclass
class DomainInteraction:
    domain: str
    substrate: str
    probability: float


def get_interaction_probabilities(interactions_file):
    domain_to_interactions = {}

    interactions_data = Tabular(interactions_file, [0, 1])
    for datapoint in interactions_data.data:
        domain_name = interactions_data.get_value(datapoint, "domain")
        substrate = interactions_data.get_value(datapoint, "substrate")
        probability = float(interactions_data.get_value(datapoint, "interaction_probability"))

        if domain_name not in domain_to_interactions:
            domain_to_interactions[domain_name] = []
        domain_to_interactions[domain_name].append(DomainInteraction(domain_name, substrate, probability))

    for domain, interactions in domain_to_interactions.items():
        interactions.sort(key=lambda x: x.probability, reverse=True)

    return domain_to_interactions


def assess_promiscuity(domain_list, interaction_file, parasect_data, included_substrates):

    promiscuous_domain_to_spec = get_promiscuous(domain_list, parasect_data, included_substrates)
    print(promiscuous_domain_to_spec)
    domain_to_interactions = get_interaction_probabilities(interaction_file)
    print(f"Analysing {len(promiscuous_domain_to_spec)} promiscuous domains..")

    correct = 0
    partially_correct = 0
    incorrect = 0

    for promiscuous_domain, spec in promiscuous_domain_to_spec.items():
        interactions = domain_to_interactions[promiscuous_domain]
        predicted_substrates = [interaction.substrate for interaction in interactions[:len(spec)]]
        print(spec, predicted_substrates)
        if set(spec) == set(predicted_substrates):
            correct += 1
        elif set(spec).intersection(set(predicted_substrates)):
            partially_correct += 1
        else:
            incorrect += 1

    print(f"{correct} correct")
    print(f"{partially_correct} partially correct")
    print(f"{incorrect} incorrect")


if __name__ == "__main__":
    # assess_promiscuity(argv[1], argv[2], argv[3], argv[4])
    plot_promiscuity(argv[1])





