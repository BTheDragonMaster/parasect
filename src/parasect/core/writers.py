from typing import Any
import os

from parasect.core.parasect_result import Result
from parasect.core.models import ModelType
from parasect.core.constants import SEPARATOR_1, SEPARATOR_2, SEPARATOR_3


def write_model_metadata_file(model_to_version: dict[ModelType, str], out_file: str) -> None:
    with open(out_file, 'w') as out:
        for model, version in model_to_version.items():
            out.write(f"sklearn_version_{model.name.lower()}\t{version}\n")


def write_results(results: list[Result],
                  out_dir: str,
                  number_predictions: int,
                  model_type: ModelType,
                  s1: str = SEPARATOR_1,
                  s2: str = SEPARATOR_2,
                  s3: str = SEPARATOR_3,
                  job_name: str = 'run_1',
                  save_signatures: bool = False,
                  save_extended_signatures: bool = False,
                  save_domain_sequences: bool = False) -> None:
    """Write PARAS results to file

    :param results: list of PARAS results
    :param out_dir: path to output directory
    :param number_predictions: number of predictions to report
    :param model_type: type of model
    :param s1: separator 1 for domain header
    :param s2: separator 2 for domain header
    :param s3: separator 3 for domain header
    :param job_name: job name
    :param save_signatures: if True, save signatures to file
    :param save_extended_signatures: if True, save extended signatures to file
    :param save_domain_sequences: if True, save domain sequences to file
    """

    if number_predictions > len(results[0].predictions):
        raise ValueError(f"Cannot report top {number_predictions}; only {len(results[0].predictions)} substrates in model")

    result_file = os.path.join(out_dir, f"{job_name}_{model_type.name.lower()}_results.txt")
    with open(result_file, 'w') as out:
        out.write("domain_id")
        for i in range(number_predictions):
            out.write(f"\tprediction_{i + 1}\tconfidence_prediction_{i + 1}")

        out.write('\n')

        for result in results:
            result.sort()
            out.write(result.get_domain_header(s1, s2, s3))
            for i in range(number_predictions):
                out.write(f"\t{result.prediction_labels[i]}\t{result.predictions[i]}")

            out.write('\n')

    id_to_sig = {}
    id_to_ext = {}
    id_to_seq = {}

    for result in results:
        domain_header = result.get_domain_header(s1, s2, s3)
        if save_signatures:
            id_to_sig[domain_header] = result.to_json()['domain_signature']
        if save_extended_signatures:
            id_to_ext[domain_header] = result.to_json()['domain_extended_signature']
        if save_domain_sequences:
            id_to_seq[domain_header] = result.to_json()['domain_sequence']

    if save_signatures:
        write_fasta_file(id_to_sig, os.path.join(out_dir, f"{job_name}_signatures.fasta"))
    if save_extended_signatures:
        write_fasta_file(id_to_ext, os.path.join(out_dir, f"{job_name}_extended_signatures.fasta"))
    if save_domain_sequences:
        write_fasta_file(id_to_seq, os.path.join(out_dir, f"{job_name}_sequences.fasta"))


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


def write_list(list_of_things: list[Any], out_file: str, sort: bool = True) -> None:
    """Write a list of things to a file, one thing per line

    :param list_of_things: list of strings
    :type list_of_things: list[Any]
    :param out_file: path to output file
    :type out_file: str
    :param sort: if True, sort list of things prior to writing
    :type sort: bool
    """
    if sort:
        list_of_things.sort()

    with open(out_file, 'w') as out:
        for thing in list_of_things:
            out.write(f"{thing}\n")
