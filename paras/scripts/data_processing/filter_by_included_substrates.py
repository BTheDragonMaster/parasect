from sys import argv

from paras.scripts.parsers.parsers import parse_domain_list, parse_specificities, parse_substrate_list


def filter_by_substrate(substrate_list, domain_list, parasect_dataset, out_file):
    substrates = parse_substrate_list(substrate_list)
    domains = parse_domain_list(domain_list)
    domain_to_specificity = parse_specificities(parasect_dataset)

    filtered_list = []
    
    for domain in domains:
        included = False
        specificities = domain_to_specificity[domain]
        for specificity in specificities:
            if specificity in substrates:
                included = True
                break

        if included:
            filtered_list.append(domain)

    with open(out_file, 'w') as out:
        for domain in filtered_list:
            out.write(f"{domain}\n")


if __name__ == "__main__":
    filter_by_substrate(argv[1], argv[2], argv[3], argv[4])

