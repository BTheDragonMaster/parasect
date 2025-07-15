from paras.scripts.parsers.parsers import parse_specificities, parse_taxonomy, parse_domain_list
from paras.scripts.parsers.tabular import Tabular
from argparse import ArgumentParser


def parse_arguments():
    parser = ArgumentParser(description="Write fungal subset of data to output file")
    parser.add_argument('-d', type=str, required=True, help="Path to file containing list of domains")
    parser.add_argument('-t', type=str, required=True, help="Path to file containing taxonomy for each domain")
    parser.add_argument('-s', type=str, required=True,
                        help="Path to file containing specificities in PARAS dataset format")
    parser.add_argument('-o', type=str, required=True, help="Path to output file")
    parser.add_argument('-u', type=str, required=True, help="Path to file containing domains of unknown taxonomy")
    parser.add_argument('-r', type=str, required=True, help="Path to file containing untrustworhty domains of unknown or fungal taxonomy")

    args = parser.parse_args()

    return args


def write_fungal_domains(domain_list_file: str, taxonomy_file: str, parasect_dataset: str, out_file: str,
                         unknowns_file: str, removed_file: str) -> None:
    domains = parse_domain_list(domain_list_file)
    taxonomy = parse_taxonomy(taxonomy_file)
    specificities = parse_specificities(parasect_dataset)
    parasect_dataset = Tabular(parasect_dataset, [0])

    with open(out_file, 'w') as out:
        with open(unknowns_file, 'w') as unknowns:
            with open(removed_file, 'w') as removed:

                out.write(f'domain_id\tsequence\tspecificity\n')
                unknowns.write('domain_id\tsequence\tspecificity\n')
                removed.write('domain_id\tsequence\tspecificity\n')
                for datapoint in parasect_dataset.data:
                    domain = parasect_dataset.get_value(datapoint, "domain_id")
                    sequence = parasect_dataset.get_value(datapoint, "sequence")
                    specificity_string = parasect_dataset.get_value(datapoint, "specificity")

                    if domain in taxonomy:

                        tax = taxonomy[domain]

                        if 'Fungi' in tax:
                            if domain in domains:
                                out.write(f"{domain}\t{sequence}\t{specificity_string}\n")

                            else:
                                removed.write(f"{domain}\t{sequence}\t{specificity_string}\n")

                    else:
                        if domain in domains:
                            unknowns.write(f"{domain}\t{sequence}\t{specificity_string}\n")
                        else:
                            removed.write(f"{domain}\t{sequence}\t{specificity_string}\n")


def main():
    args = parse_arguments()
    write_fungal_domains(args.d, args.t, args.s, args.o, args.u, args.r)


if __name__ == "__main__":
    main()
