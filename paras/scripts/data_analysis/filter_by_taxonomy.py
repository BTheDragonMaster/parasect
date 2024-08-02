from paras.scripts.parsers.parsers import parse_specificities, parse_taxonomy, parse_domain_list
from paras.scripts.parsers.tabular import Tabular
from sys import argv


taxonomy = parse_taxonomy(argv[1])
dataset = Tabular(argv[2], [0])
out_file = argv[3]

with open(out_file, 'w') as out:
    out.write("domain_id\tsequence\tspecificity\n")
    for domain in dataset.data:
        domain_id = dataset.get_value(domain, "domain_id")
        sequence = dataset.get_value(domain, "sequence")
        specificity = dataset.get_value(domain, "specificity")
        protein = '.'.join(domain_id.split('.')[:-1])

        if protein in taxonomy:
            tax = taxonomy[protein]
            if 'Fungi' in tax:
                out.write(f"{domain_id}\t{sequence}\t{specificity}\n")