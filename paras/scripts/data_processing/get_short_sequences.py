from argparse import ArgumentParser, Namespace
from paras.scripts.parsers.tabular import Tabular
from paras.scripts.parsers.parsers import parse_domain_list


def parse_args() -> Namespace:
    parser = ArgumentParser(description="Get sequences of specified length")
    parser.add_argument("-i", required=True, type=str, help="Path to parasect dataset")
    parser.add_argument("-o", required=True, type=str, help="Path to output file")
    parser.add_argument("-ut", required=True, type=int, help="Upper threshold (inclusive)")
    parser.add_argument("-lt", type=int, default=0, help="Lower threshold (inclusive)")
    parser.add_argument("-d", type=str, required=True, help="Path to file containing domain list")

    args = parser.parse_args()
    return args


def get_sequences_by_length(parasect_dataset: str, domain_list: str, out_file: str,
                            upper_threshold: int, lower_threshold: int) -> None:
    parasect_data = Tabular(parasect_dataset, [0])
    domains = parse_domain_list(domain_list)

    with open(out_file, 'w') as out:
        out.write(f"domain_id\tsequence\tspecificity\n")
        for datapoint in parasect_data.data:
            domain_id = parasect_data.get_value(datapoint, "domain_id")
            if domain_id in domains:
                sequence = parasect_data.get_value(datapoint, "sequence")
                specificity_string = parasect_data.get_value(datapoint, "specificity")
                if lower_threshold <= len(sequence) <= upper_threshold:
                    out.write(f"{domain_id}\t{sequence}\t{specificity_string}\n")


def main():
    args = parse_args()
    get_sequences_by_length(args.i, args.d, args.o, args.ut, args.lt)


if __name__ == "__main__":
    main()
