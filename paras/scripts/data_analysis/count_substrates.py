from collections import OrderedDict
from paras.scripts.parsers.parsers import parse_specificities, parse_domain_list, parse_substrate_list


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
