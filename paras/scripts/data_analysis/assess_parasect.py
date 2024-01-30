from sys import argv

from paras.scripts.parsers.tabular import Tabular
import matplotlib.pyplot as plt


def count_correct_in_top_x(interaction_file, top_x):
    interaction_data = Tabular(interaction_file, [0, 1])
    domain_to_interactions = {}

    for datapoint in interaction_data.data:
        domain = interaction_data.get_value(datapoint, "domain")
        probability = interaction_data.get_value(datapoint, "interaction_probability")
        true_interaction = interaction_data.get_value(datapoint, "true_interaction")
        if true_interaction == 'True':
            true_interaction = True
        else:
            true_interaction = False

        if domain not in domain_to_interactions:
            domain_to_interactions[domain] = []

        domain_to_interactions[domain].append((probability, true_interaction))

    total_correct = 0
    total_incorrect = 0

    for domain, interactions in domain_to_interactions.items():
        interactions.sort(reverse=True)
        top_interactions = interactions[:top_x]
        correct = False
        for probability, true_interaction in top_interactions:
            if true_interaction:
                correct = True

        if correct:
            total_correct += 1
        else:
            total_incorrect += 1

    accuracy = total_correct / (total_correct + total_incorrect)

    print(f"Frequency of substrates occurring in {top_x}:")

    print(f"Correct: {total_correct}")
    print(f"Incorrect: {total_incorrect}")
    print(f"Accuracy: {accuracy}")

    return accuracy


def plot_interactions(interaction_file, top_range, out_file):
    x = []
    y = []

    for i in range(top_range):
        accuracy = count_correct_in_top_x(interaction_file, i + 1)
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
    plt.xticks([1, 2, 3, 4, 5], [1, 2, 3, 4, 5])
    plt.xlim(1, max(x) + 1)
    plt.ylim(70, 100)

    plt.savefig(out_file, bbox_inches='tight')


if __name__ == "__main__":
    plot_interactions(argv[1], int(argv[2]), argv[3])

