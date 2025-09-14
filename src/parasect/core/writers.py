from parasect.core.parasect_result import Result


def write_results(results: list[Result], out_file: str, number_predictions: int) -> None:
    """Write PARAS results to file

    :param results: list of PARAS results
    :type results: list[Result]
    :param out_file: path to output file
    :type out_file: str
    :param number_predictions: number of predictions to report
    :type number_predictions: int
    """

    if number_predictions > len(results[0].predictions):
        raise ValueError(f"Cannot report top {number_predictions}; only {len(results[0].predictions)} substrates in model")

    with open(out_file, 'w') as out:
        out.write("domain_id")
        for i in range(number_predictions):
            out.write(f"\tprediction_{i + 1}\tconfidence_prediction_{i + 1}")

        out.write('\n')

        for result in results:
            result.sort()
            out.write(result.get_domain_header())
            for i in range(number_predictions):
                out.write(f"\t{result.prediction_labels[i]}\t{result.predictions[i]}")

            out.write('\n')


def write_fasta_file(fasta_dict: dict[str, str], path_out: str) -> None:
    """Write a dictionary of fasta sequences to a file.

    :param fasta_dict: Dictionary of fasta sequences, where the key is the sequence
        header and the value is the sequence.
    :type fasta_dict: Dict[str, str]
    :param path_out: Path to output fasta file.
    :type path_out: str
    """
    sorted_ids = sorted(fasta_dict.keys())
    with open(path_out, "w") as fo:

        # iterate over the dictionary items
        for header in sorted_ids:
            sequence = fasta_dict[header]
            fo.write(f">{header}\n{sequence}\n")


def write_list(list_of_things: list[str], out_file: str) -> None:
    """Write a list of things to a file, one thing per line

    :param list_of_things: list of strings
    :type list_of_things: list[str]
    :param out_file: path to output file
    :type out_file: str
    """

    list_of_things.sort()

    with open(out_file, 'w') as out:
        for thing in list_of_things:
            out.write(f"{thing}\n")
