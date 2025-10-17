from argparse import ArgumentParser, Namespace
from typing import Optional


from parasect.core.parsing import parse_raw_taxonomy
from parasect.core.tabular import Tabular


def parse_arguments() -> Namespace:
    """Parse arguments from command line

    :return: Arguments
    :rtype: Namespace
    """
    parser = ArgumentParser(description="Retrieve taxonomy for database domains")
    parser.add_argument("-r", '--rank_file', required=True, type=str,
                        help="Path to PARASECT database")
    parser.add_argument("-o", '--out', required=True, type=str,
                        help="Path to output file")
    parser.add_argument("-t", '--taxonomy_file', required=True, type=str,
                        help="Path to taxonomy file")
    parser.add_argument("-m", "--maintain", nargs="+", type=str, help="Ranks to maintain")

    arguments = parser.parse_args()

    return arguments


def get_condensed_taxonomy(taxonomy_file: str, rank_mapping: str, ranks_to_maintain: list[str], out_file: str) -> None:
    """Condense taxonomy to only contain predefined ranks

    :param taxonomy_file: path to taxonomy file
    :type taxonomy_file: str
    :param rank_mapping: path to file containing a mapping from taxonomic name to rank
    :type rank_mapping: str
    :param ranks_to_maintain: ranks to maintain in condensed taxonomy
    :type ranks_to_maintain: list[str]
    :param out_file: path to output file
    :type out_file: str
    """

    rank_data = Tabular(rank_mapping)

    id_to_tax = parse_raw_taxonomy(taxonomy_file)
    id_to_condensed_tax: dict[str, dict[str, Optional[str]]] = {}

    for prot_id, tax in id_to_tax.items():
        id_to_condensed_tax[prot_id] = {}
        for rank in ranks_to_maintain:
            id_to_condensed_tax[prot_id][rank] = None

        for tax_element in tax:
            rank = rank_data.get_row_value(tax_element, "rank")
            if rank in ranks_to_maintain:
                id_to_condensed_tax[prot_id][rank] = tax_element

    with open(out_file, 'w') as out:
        out.write('protein_id')
        for rank in ranks_to_maintain:
            out.write(f"\t{rank}")
        out.write('\n')
        proteins = sorted(id_to_condensed_tax.keys())
        for protein in proteins:
            out.write(f"{protein}")
            for rank in ranks_to_maintain:
                out.write(f"\t{id_to_condensed_tax[protein][rank]}")
            out.write('\n')


def main():
    args = parse_arguments()
    get_condensed_taxonomy(args.taxonomy_file, args.rank_file, args.maintain, args.out)


if __name__ == "__main__":
    main()
