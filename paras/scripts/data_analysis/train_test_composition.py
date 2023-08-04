from sys import argv

from paras.scripts.parsers.parsers import parse_specificities, parse_domain_list
from paras.scripts.data_analysis.count_substrates import count_substrates


def assess_train_test_composition(train, test, specificities, included_specs_file):
    included_substrates = []
    with open(included_specs_file, 'r') as included_specs:
        for line in included_specs:
            included_substrates.append(line.strip())

    domain_list_train = parse_domain_list(train)
    domain_list_test = parse_domain_list(test)

    train_specs = []
    test_specs = []

    domain_to_spec = parse_specificities(specificities)

    for domain in domain_list_train:
        train_specs += domain_to_spec[domain]

    train_specs.sort()

    for domain in domain_list_test:
        test_specs += domain_to_spec[domain]

    test_specs.sort()

    train_counts = count_substrates(train_specs)
    test_counts = count_substrates(test_specs)

    for substrate, train_count in train_counts.items():
        if substrate in included_substrates:
            test_count = test_counts[substrate]
            print(f"{substrate}: {test_count / (train_count + test_count)}")


if __name__ == "__main__":
    assess_train_test_composition(argv[1], argv[2], argv[3], argv[4])
