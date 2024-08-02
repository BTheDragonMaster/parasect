from paras.scripts.parsers.parsers import parse_specificities, parse_taxonomy, parse_domain_list
from sys import argv

domains = parse_domain_list(argv[1])
taxonomy = parse_taxonomy(argv[2])
specificities = parse_specificities(argv[3])

for domain in domains:
    print(domain)
    protein = '.'.join(domain.split('.')[:-1])
    if protein in taxonomy:
        print('\t'.join(taxonomy[protein][:3]))
    print(specificities[domain])
    print('\n')