from paras.scripts.parsers.parsers import parse_specificities, parse_taxonomy, parse_domain_list
from sys import argv

domains = parse_domain_list(argv[1])
taxonomy = parse_taxonomy(argv[2])
specificities = parse_specificities(argv[3])

counter = 0
for domain in domains:

    protein = '.'.join(domain.split('.')[:-1])
    if protein in taxonomy:
        tax = taxonomy[protein]
        if 'Fungi' in tax:
            counter += 1
            print(domain)
            print('\t'.join(specificities[domain]))
            print('\t'.join(taxonomy[protein][:3]))


print(counter)