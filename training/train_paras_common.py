# -*- coding: utf-8 -*-

"""Script for training the PARAS model for predicting common substrate specificities."""

import argparse 
import logging


def cli() -> argparse.Namespace:
    """Command line interface for training PARAS models."""
    parser = argparse.ArgumentParser()

    # required arguments
    parser.add_argument(
        "--data-dir",
        type=str,
        required=True,
        help="path to the data directory",
    )

    parser.add_argument(
        "--out-dir",
        type=str,
        required=True,
        help="path to the output directory",
    )

    # optional arguments
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        required=False,
        default=1,
        help="number of threads to use",
    )

    parser.add_argument(
        "-l",
        "--logger-level",
        type=str,
        required=False,
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="set the logging level",
    )

    return parser.parse_args()


def main() -> None:
    """Run the training script."""
    args = cli()
    
    logging.basicConfig(level=args.logger_level)
    logger = logging.getLogger(__name__)

    logger.debug(f"command line arguments: {args}")

    logger.error("Not implemented yet.")
    exit(1)


if __name__ == "__main__":
    main()
