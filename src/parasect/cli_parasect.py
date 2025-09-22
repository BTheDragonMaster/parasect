# -*- coding: utf-8 -*-

"""CLI for PARASECT."""

import os
import argparse
import logging
from joblib import load

from parasect.core.constants import SEPARATOR_1, SEPARATOR_2, SEPARATOR_3
from parasect.api import run_parasect
from parasect.core.helpers import download_and_unpack_or_fetch
from parasect.core.writers import write_fasta_file, write_results
from parasect.core.parsing import parse_smiles_mapping


def cli() -> argparse.Namespace:
    """CLI for PARASECT.

    :return: CLI arguments
    :rtype: argparse.Namespace
    """

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', "--input", type=str, required=True, help="Path to input fasta or gbk file.")
    parser.add_argument('-f', "--file_type", type=str, default='fasta',
                        help="Input file type. Must be 'fasta' or 'gbk'.")
    parser.add_argument('-o', "--output", type=str, required=True, help="Path to output directory.")
    parser.add_argument('-j', "--job_name", type=str, default="run_1",
                        help="Job name")
    parser.add_argument('-n', "--number_predictions", type=int, default=3, help="Number of top predictions to report.")
    parser.add_argument('-t', "--temp", type=str, default=None,
                        help="Temp dir. If not given, create temp folder in output dir")
    parser.add_argument('-p', "--profile_alignment", action='store_true',
                        help="Use profile alignment instead of HMM for active site extraction")
    parser.add_argument('-save_extended', action='store_true',
                        help="Save extended 34 amino acid signatures to file.")
    parser.add_argument('-save_signatures', action='store_true',
                        help="Save short 10 amino acid signatures to file ('Stachelhaus code')")
    parser.add_argument('-save_domains', action='store_true',
                        help="Save full a domain sequences to file ('Stachelhaus code')")

    parser.add_argument('-s1', type=str, default=SEPARATOR_1, help="Symbol used as separator")
    parser.add_argument('-s2', type=str, default=SEPARATOR_2, help="Symbol used as separator")
    parser.add_argument('-s3', type=str, default=SEPARATOR_3, help="Symbol used as separator")

    parser.add_argument('-exclude_standard_substrates', action='store_true',
                        help="Don't run predictions for the default substrates included in PARASECT")
    parser.add_argument('-smiles', type=str, default=None,
                        help="File containing custom substrate names (column 1) and SMILES (column 2), with a header")

    args = parser.parse_args()
    assert args.file_type.upper() in ["FASTA", "GBK"]

    return args


def main() -> None:
    """Run CLI for PARASECT."""
    """Run CLI for PARAS."""
    args = cli()
    logger = logging.getLogger(__name__)
    logging.basicConfig(level="INFO")

    if not os.path.exists(args.output):
        os.mkdir(args.output)

    if args.temp is None:
        temp_dir = os.path.join(args.output, "temp")
        if not os.path.exists(temp_dir):
            os.mkdir(temp_dir)

    else:
        temp_dir = args.temp

    if args.smiles is not None:
        substrates = parse_smiles_mapping(args.smiles)
        substrate_names = [s.name for s in substrates]
        substrate_smiles = [s.smiles for s in substrates]
    else:
        substrate_names = None
        substrate_smiles = None

    if not substrate_names and args.exclude_standard_substrates:
        raise ValueError("No substrates to test! Either include standard substrates or pass custom substrate SMILES")

    model_path = download_and_unpack_or_fetch(
        r"https://zenodo.org/records/17155186/files/model.parasect.gz?download=1",
        temp_dir, logger)

    model = load(model_path)

    with open(args.input, 'r') as input_file:
        protein_data = input_file.read()

    results = run_parasect(protein_data, args.file_type, temp_dir, model,
                           custom_substrate_names=substrate_names,
                           custom_substrate_smiles=substrate_smiles,
                           only_custom=args.exclude_standard_substrates,
                           use_structure_guided_alignment=args.profile_alignment)

    id_to_sig = {}
    id_to_ext = {}
    id_to_seq = {}

    for result in results:
        domain_header = result.get_domain_header(args.s1, args.s2, args.s3)
        if args.save_signatures:
            id_to_sig[domain_header] = result.to_json()['domain_signature']
        if args.save_extended:
            id_to_ext[domain_header] = result.to_json()['domain_extended_signature']
        if args.save_domains:
            id_to_seq[domain_header] = result.to_json()['domain_sequence']

    if args.save_signatures:
        write_fasta_file(id_to_sig, os.path.join(args.output, f"{args.job_name}_signatures.fasta"))
    if args.save_extended:
        write_fasta_file(id_to_ext, os.path.join(args.output, f"{args.job_name}_extended_signatures.fasta"))
    if args.save_domains:
        write_fasta_file(id_to_seq, os.path.join(args.output, f"{args.job_name}_sequences.fasta"))

    results_out = os.path.join(args.output, f"{args.job_name}_parasect_results.txt")
    write_results(results, results_out, args.number_predictions)


if __name__ == "__main__":
    main()
