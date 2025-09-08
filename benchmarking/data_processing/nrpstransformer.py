from argparse import ArgumentParser, Namespace

from parasect.core.tabular import Tabular, write_tabular


def parse_arguments() -> Namespace:
    """
    Parse command line arguments

    :return: arguments
    :rtype: Namespace
    """
    parser = ArgumentParser(description="Process nrpstransformer data into parasect dataset-type file")
    parser.add_argument('-i', type=str, required=True,
                        help="Path to input file containing raw nrpstransformer data in tab-separated format")
    parser.add_argument('-m', type=str, required=True,
                        help="Path to file containing substrate mapping from NRPSTransformer abbreviation to parasect substrate(s)")
    parser.add_argument('-o', type=str, required=True,
                        help="Path to output file")

    args = parser.parse_args()

    return args


def process_nrpstransformer_data(nrpstransformer_raw_file: str, substrate_mapping_file: str, out_file: str) -> None:
    """
    Write nrpstransformer data in parasect format

    :param nrpstransformer_raw_file: file containing raw nrpstransformer data
    :type nrpstransformer_raw_file: str
    :param substrate_mapping_file: file containing nrpstransformer substrate mapping
    :type substrate_mapping_file: str
    :param out_file: output file
    :type out_file: str
    """

    id_to_seq: dict[str, str] = {}
    id_to_spec: dict[str, str] = {}
    abbreviation_to_spec: dict[str, str] = {}

    substrate_mapping = Tabular(substrate_mapping_file)
    for abbreviation in substrate_mapping.rows:
        specificity = substrate_mapping.get_row_value(abbreviation, "full_name")
        abbreviation_to_spec[abbreviation] = specificity

    raw_data = Tabular(nrpstransformer_raw_file)
    for domain_name in raw_data.rows:
        if "sequence for NRPStransformer development" in raw_data.column_names:
            sequence = raw_data.get_row_value(domain_name, "sequence for NRPStransformer development")
        elif "A domain Sequence full length" in raw_data.column_names:
            sequence = raw_data.get_row_value(domain_name, "A domain Sequence full length")
        else:
            raise ValueError(f"No sequence column found in data file")

        if "Label amino acid" in raw_data.column_names:
            specificity = raw_data.get_row_value(domain_name, "Label amino acid")
        elif "Label" in raw_data.column_names:
            specificity = raw_data.get_row_value(domain_name, "Label")
        else:
            raise ValueError(f"No label column found in data file")

        id_to_seq[domain_name] = sequence
        specificities = specificity.lower().split('|')

        specificity_string = []
        for specificity in specificities:
            spec = abbreviation_to_spec[specificity.lower()]
            specs = spec.split('|')
            specificity_string.extend(specs)
        specificity_string = list(set(specificity_string))
        id_to_spec[domain_name] = '|'.join(specificity_string)

    write_tabular([id_to_seq, id_to_spec], ["domain_id", "sequence", "specificity"], out_file)


if __name__ == "__main__":
    args = parse_arguments()
    process_nrpstransformer_data(args.i, args.m, args.o)
