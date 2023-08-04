from paras.scripts.parsers.parsers import parse_substrate_list


def relabel_from_inclusion_limit(domain_list, domain_to_substrates, substrate_counts, limit):
    domain_to_filtered = {}
    for domain in domain_list:
        for substrate in domain_to_substrates[domain]:
            if substrate_counts[substrate] >= limit:
                if domain not in domain_to_filtered:
                    domain_to_filtered[domain] = []
                domain_to_filtered[domain].append(substrate)

    return domain_to_filtered


def relabel_from_substrate_list(domain_list, domain_to_substrates, substrate_file):
    included_substrates = parse_substrate_list(substrate_file)

    domain_to_filtered = {}
    for domain in domain_list:
        for substrate in domain_to_substrates[domain]:
            if substrate in included_substrates:
                if domain not in domain_to_filtered:
                    domain_to_filtered[domain] = []
                domain_to_filtered[domain].append(substrate)

    return domain_to_filtered
