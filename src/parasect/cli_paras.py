# -*- coding: utf-8 -*-

"""CLI for PARAS."""

import os
import argparse
import logging
from joblib import load
from shutil import copy

from parasect.core.constants import SEPARATOR_1, SEPARATOR_2, SEPARATOR_3
from parasect.api import run_paras
from parasect.core.helpers import download_and_unpack_or_fetch
from parasect.core.writers import write_fasta_file, write_results
from parasect.core.retrain_models import retrain_model, model_needs_retraining, update_metadata_file
from parasect.core.models import ModelType
from parasect.core.constants import MODEL_METADATA_FILE


def cli() -> argparse.Namespace:
    """CLI for PARAS.

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
    parser.add_argument('-t', "--temp", type=str, default=None, help="Temp dir. If not given, create temp folder in output dir")
    parser.add_argument('-p', "--profile_alignment", action='store_true',
                        help="Use profile alignment instead of HMM for active site extraction")
    parser.add_argument('-save_extended', action='store_true',
                        help="Save extended 34 amino acid signatures to file.")
    parser.add_argument('-save_signatures', action='store_true',
                        help="Save short 10 amino acid signatures to file ('Stachelhaus code')")
    parser.add_argument('-save_domains', action='store_true',
                        help="Save full a domain sequences to file ('Stachelhaus code')")
    parser.add_argument('-all_substrates', action='store_true',
                        help="Use classifier that was trained on all substrates")
    parser.add_argument('-s1', type=str, default=SEPARATOR_1, help="Symbol used as separator")
    parser.add_argument('-s2', type=str, default=SEPARATOR_2, help="Symbol used as separator")
    parser.add_argument('-s3', type=str, default=SEPARATOR_3, help="Symbol used as separator")

    args = parser.parse_args()
    assert args.file_type.upper() in ["FASTA", "GBK"]

    return args


def main() -> None:
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

    metadata_path = os.path.join(temp_dir, "model_metadata.txt")
    if not os.path.exists(metadata_path):
        copy(MODEL_METADATA_FILE, metadata_path)

    if args.all_substrates:
        if model_needs_retraining(metadata_path, ModelType.PARAS_ALL_SUBSTRATES):
            print("Found incompatible version of scikit-learn. Retraining..")
            model = retrain_model(ModelType.PARAS_ALL_SUBSTRATES)
            model_path = os.path.join(temp_dir, model.file_name)
            model.save(temp_dir)
            update_metadata_file(ModelType.PARAS_ALL_SUBSTRATES, metadata_path)

        else:

            model_path = download_and_unpack_or_fetch(r"https://zenodo.org/records/17224548/files/all_substrates_model.paras.gz?download=1",
                                                      temp_dir, logger)

    else:

        if model_needs_retraining(metadata_path, ModelType.PARAS):
            print("Found incompatible version of scikit-learn. Retraining..")
            model = retrain_model(ModelType.PARAS)
            model_path = os.path.join(temp_dir, model.file_name)
            model.save(temp_dir)
            update_metadata_file(ModelType.PARAS, metadata_path)

        else:
            model_path = download_and_unpack_or_fetch(
                r"https://zenodo.org/records/17224548/files/model.paras.gz?download=1",
                temp_dir, logger)

    model = load(model_path)

    with open(args.input, 'r') as input_file:
        protein_data = input_file.read()

    results = run_paras(protein_data, args.file_type, temp_dir, model, args.profile_alignment)

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

    results_out = os.path.join(args.output, f"{args.job_name}_paras_results.txt")
    write_results(results, results_out, args.number_predictions)


if __name__ == "__main__":
    main()
