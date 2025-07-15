from collections import OrderedDict
from paras.scripts.parsers.parsers import parse_specificities, parse_domain_list, parse_substrate_list, parse_substrate_smiles
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

def write_counts(domain_file, specificity_file, output_file):
    spec_to_count = count_substrates_from_file(domain_file, specificity_file)
    with open(output_file, 'w') as out:
        out.write("substrate\tcount\n")
        specs = list(spec_to_count.keys())
        specs.sort(key = lambda x: spec_to_count[x], reverse=True)
        for spec in specs:
            out.write(f"{spec}\t{spec_to_count[spec]}\n")

def write_cytoscape_sizes(domain_file, specificity_file, smiles_file, output_file, out_smiles):
    spec_to_smiles = parse_substrate_smiles(smiles_file)
    spec_to_count = count_substrates_from_file(domain_file, specificity_file)
    with open(output_file, 'w') as out:
        with open(out_smiles, 'w') as out_2:
            out_2.write(f"substrate\tsmiles\n")
            specs = list(spec_to_count.keys())
            specs.sort(key=lambda x: spec_to_count[x], reverse=True)
            for spec, count in spec_to_count.items():
                out_2.write(f"{spec}\t{spec_to_smiles[spec]}\n")
                if count > 100:
                    out.write(f"{spec}\t4\n")
                elif count > 50:
                    out.write(f"{spec}\t3\n")
                elif count > 10:
                    out.write(f"{spec}\t2\n")
                else:
                    out.write(f"{spec}\t1\n")

    print('\n')
    for spec, smiles in spec_to_smiles.items():
        if spec not in spec_to_count:
            print(spec)




if __name__ == "__main__":
    # get_included_substrates(argv[1], argv[2], argv[3])
    write_cytoscape_sizes(argv[1], argv[2], argv[4], argv[3], argv[5])