# -*- coding: utf-8 -*-

"""Example inputs users can load on the submit and annotation pages."""

from __future__ import annotations

import os

from flask import Blueprint, Response

from .common import ResponseData, Status

EXAMPLE_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "example_data")

# Only files listed here can be served, so the ID in the URL never touches the
# file system. To add an example:
#       - drop the file into example_data/
#       - add an entry here
EXAMPLE_INPUTS: dict[str, dict[str, str]] = {
    "dptA": {
        "label": "dptA protein (FASTA)",
        "description": "Daptomycin NRPS dptA from Streptomyces roseosporus, a single protein with several A domains.",
        "fileName": "dptA.fasta",
        "inputType": "fasta",
    },
    "BGC0000336": {
        "label": "Daptomycin gene cluster (GenBank)",
        "description": "MIBiG BGC0000336 from Streptomyces filamentosus NRRL 11379, a full biosynthetic gene cluster.",
        "fileName": "BGC0000336.gbk",
        "inputType": "gbk",
    },
}


def read_example_input(example_id: str) -> str:
    """Read the contents of an example input file.

    :param example_id: Key in EXAMPLE_INPUTS.
    :return: File contents.
    :raises KeyError: If the example does not exist.
    """
    with open(os.path.join(EXAMPLE_DATA_DIR, EXAMPLE_INPUTS[example_id]["fileName"])) as f:
        return f.read()


blueprint_example_inputs = Blueprint("example_inputs", __name__)


@blueprint_example_inputs.route("/api/example_inputs", methods=["GET"])
def list_example_inputs() -> Response:
    """List the available example inputs, without their contents."""
    examples = [{"id": example_id, **meta} for example_id, meta in EXAMPLE_INPUTS.items()]
    return ResponseData(Status.Success, payload={"examples": examples}).to_dict()


@blueprint_example_inputs.route("/api/example_inputs/<example_id>", methods=["GET"])
def get_example_input(example_id: str) -> Response:
    """Return one example input, including its contents."""
    if example_id not in EXAMPLE_INPUTS:
        return ResponseData(Status.Failure, message=f"unknown example: {example_id}").to_dict()
    try:
        content = read_example_input(example_id)
    except OSError as e:
        return ResponseData(Status.Failure, message=f"failed to read example: {str(e)}").to_dict()
    return ResponseData(
        Status.Success,
        payload={"id": example_id, **EXAMPLE_INPUTS[example_id], "content": content},
    ).to_dict()
