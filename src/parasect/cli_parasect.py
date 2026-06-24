# -*- coding: utf-8 -*-

"""CLI for PARASECT."""

import os
import argparse
import logging
from joblib import load
from shutil import rmtree

from parasect.core.constants import SEPARATOR_1, SEPARATOR_2, SEPARATOR_3
from parasect.api import run_parasect
from parasect.core.helpers import prepare_folders, prepare_substrates, prepare_model
from parasect.core.writers import write_fasta_file, write_results
from parasect.core.models import ModelType


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
    parser.add_argument('-n', "--number_predictions", type=int, default=3,
                        help="Number of top predictions to report.")
    parser.add_argument('-t', "--temp", type=str, default=None,
                        help="Temp dir. If not given, create temp folder in output dir")
    parser.add_argument('-p', "--profile_alignment", action='store_true',
                        help="Use profile alignment instead of HMM for active site extraction")
    parser.add_argument('-m', "--model_dir", type=str, default=None,
                        help="Path to model directory. If not given, use temp folder")
    parser.add_argument('-save_extended', action='store_true',
                        help="Save extended 34 amino acid signatures to file.")
    parser.add_argument('-save_signatures', action='store_true',
                        help="Save short 10 amino acid signatures to file ('Stachelhaus code')")
    parser.add_argument('-save_domains', action='store_true',
                        help="Save full a domain sequences to file ('Stachelhaus code')")
    parser.add_argument('-bacterial', action='store_true',
                        help="If given, run bacterial-only model")

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
    args = cli()
    logger = logging.getLogger(__name__)
    logging.basicConfig(level="INFO")

    temp_dir, model_dir = prepare_folders(args.output, args.temp, args.model_dir)

    if args.bacterial:
        model_type = ModelType.PARASECT_BACTERIAL
    else:
        model_type = ModelType.PARASECT

    substrate_names, substrate_smiles = prepare_substrates(args.smiles)

    if not substrate_names and args.exclude_standard_substrates:
        raise ValueError("No substrates to test! Either include standard substrates or pass custom substrate SMILES")

    model_path = prepare_model(model_type, model_dir, logger)
    model = load(model_path)

    with open(args.input, 'r') as input_file:
        protein_data = input_file.read()

    results = run_parasect(protein_data, args.file_type, temp_dir, model,
                           custom_substrate_names=substrate_names,
                           custom_substrate_smiles=substrate_smiles,
                           only_custom=args.exclude_standard_substrates,
                           use_structure_guided_alignment=args.profile_alignment,
                           bacterial_only=args.bacterial)

    write_results(results, args.output, args.number_predictions, model_type,
                  args.s1, args.s2, args.s3,
                  args.job_name,
                  args.save_signatures,
                  args.save_extended,
                  args.save_domains)

    rmtree(temp_dir)


if __name__ == "__main__":
    main()
