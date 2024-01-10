from collections import OrderedDict
from paras.scripts.parsers.parsers import parse_specificities, parse_domain_list, parse_substrate_list
from sys import argv

def count_substrates(spec_list):
    counts = OrderedDict()
    for specificity in spec_list:
        if specificity not in counts:
            counts[specificity] = 0
        counts[specificity] += 1

    return counts


def count_substrates_from_file(domain_file, specificity_file, substrate_file=None):
    included_substrates = []
    if substrate_file is not None:
        included_substrates = parse_substrate_list(substrate_file)

    domain_list = parse_domain_list(domain_file)
    domain_to_spec = parse_specificities(specificity_file)

    all_specs = []
    for domain in domain_list:
        all_specs += domain_to_spec[domain]

    counts = count_substrates(all_specs)
    if substrate_file is not None:
        for substrate in list(counts.keys()):
            if substrate not in included_substrates:
                del counts[substrate]
    return counts


def get_included_substrates(domain_file, specificity_file, output_file, min_count=1):
    domain_list = parse_domain_list(domain_file)
    domain_to_spec = parse_specificities(specificity_file)

    all_specs = []
    for domain in domain_list:
        all_specs += domain_to_spec[domain]

    counts = count_substrates(all_specs)
    with open(output_file, 'w') as out:
        for substrate, count in counts.items():
            if count >= min_count:
                out.write(f"{substrate}\n")

if __name__ == "__main__":
    get_included_substrates(argv[1], argv[2], argv[3])