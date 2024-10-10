# -*- coding: utf-8 -*-

"""Module for handling MUSCLE requests."""

import subprocess


def run_muscle(path_in: str, path_in_alignment: str, path_out: str) -> None:
    """Run MUSCLE with the given input file and alignment file.

    :param path_in: path to the input file
    :type path_in: str
    :param path_in_alignment: path to the alignment file
    :type path_in_alignment: str
    :param path_out: path to the output file
    :type path_out: str
    """
    # compile the command for running MUSCLE
    command = [
        "muscle",
        "-quiet",
        "-profile",
        "-in1",
        path_in_alignment,
        "-in2",
        path_in,
        "-out",
        path_out,
    ]

    # use subprocess to run the compiled command
    subprocess.check_call(
        command,
        stdout=subprocess.DEVNULL,  # hide stdout
        stderr=subprocess.DEVNULL,  # hide stderr
    )
