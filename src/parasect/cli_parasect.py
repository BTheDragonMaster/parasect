# -*- coding: utf-8 -*-

"""CLI for PARASECT."""

import argparse
import logging


def cli() -> argparse.Namespace:
    """CLI for PARASECT.

    :return: CLI arguments
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-l",
        "--logger-level",
        required=False,
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="set the logger level",
    )

    return parser.parse_args()


def main() -> None:
    """Run CLI for PARASECT."""
    args = cli()

    logging.basicConfig(level=args.logger_level)
    logger = logging.getLogger(__name__)

    logger.error("command line interface not available for the webapp branch")


if __name__ == "__main__":
    main()
