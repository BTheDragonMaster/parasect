import os
from argparse import ArgumentParser, Namespace

from parasect.core.tabular import Tabular


def parse_arguments() -> Namespace:
    """Parse arguments from command line

    :return: Arguments
    :rtype: Namespace
    """
    parser = ArgumentParser(description="Retrieve taxonomy for database domains")
    parser.add_argument("-i", '--input', required=True, type=str,
                        help="Path to file listing substrates for each taxonomic rank")
    parser.add_argument("-s", '--substrate_summary', required=True, type=str,
                        help="Path to file listing substrates for each taxonomic rank")
    parser.add_argument("-c", '--cutoff', default=6, type=int,
                        help="Minimum substrate count for inclusion")

    arguments = parser.parse_args()

    return arguments


def parse_counts(count_file: str, cutoff: int = 1) -> dict[str, int]:
    substrate_to_count: dict[str, int] = {}
    with open(count_file, 'r') as counts:
        for line in counts:
            line = line.strip()
            if line:
                substrate, count = line.split('\t')
                count = int(count)
                if count >= cutoff:
                    substrate_to_count[substrate] = count

    return substrate_to_count


def explore_data_split(substrate_count_file: str, rank_count_file: str, cutoff: int) -> None:
    substrate_to_count = parse_counts(substrate_count_file, cutoff)
    print(f"Number of included substrates: {len(substrate_to_count)}")
    rank_data = Tabular(rank_count_file)
    for substrate in rank_data.rows:
        if substrate in substrate_to_count:
            number_of_clades_with_substrate = 0
            for category in rank_data.column_names[1:]:
                count = int(rank_data.get_row_value(substrate, category))
                if count > 0:
                    number_of_clades_with_substrate += 1

            if number_of_clades_with_substrate < 2:
                print(f"Splitting on this rank not possible for substrate: {substrate}")


def main():
    args = parse_arguments()
    explore_data_split(args.substrate_summary, args.input, args.cutoff)


if __name__ == "__main__":
    main()
