# -*- coding: utf-8 -*-

"""Routes for making adenylation domain subtrate specificity predictions with PARASECT."""

import os

import joblib
from flask import Blueprint, Response, request

from parasect.api import run_parasect
from parasect.core.constants import FINGERPRINTS_FILE, INCLUDED_SUBSTRATES_FILE
from parasect.core.parsing import (
    bitvector_from_smiles,
    bitvectors_from_substrate_names,
    parse_substrate_list,
)

from .common import ResponseData, Status
from .constants import MODEL_DIR, TEMP_DIR

blueprint_submit_parasect = Blueprint("submit_parasect", __name__)


@blueprint_submit_parasect.route("/api/submit_parasect", methods=["POST"])
def submit_parasect() -> Response:
    """Submit settings for prediction with Parasect model.

    :return: Response.
    :rtype: Response
    """
    data = request.get_json()

    # read settings
    # is everything present? return with message if not
    try:
        data = data["data"]
        selected_input = data["src"]  # fasta or gbk file contents
        selected_input_type = data["selectedInputType"]  # fasta or gbk
        num_predictions_to_report = data["numPredictionsToReport"]
        smiles_input = data["smilesSrc"]
        only_predict_for_smiles_input = data["onlyMakePredictionsUploadedSmiles"]
        use_structure_guided_alignment = data["useStructureGuidedProfileAlignment"]
        first_separator = data["firstSeparator"]
        second_separator = data["secondSeparator"]
        third_separator = data["thirdSeparator"]

    except Exception as e:
        msg = f"failed to read settings: {str(e)}"
        return ResponseData(Status.Failure, message=msg).to_dict()

    # sanity check selected input type
    # return with message if invalid
    selected_input_type = selected_input_type.strip().lower()
    if selected_input_type not in ["fasta", "gbk"]:
        msg = f"invalid input type: {selected_input_type}."
        return ResponseData(Status.Failure, message=msg).to_dict()

    # locate temp directory
    # return error if not found.
    try:
        assert os.path.exists(TEMP_DIR)
    except Exception as e:
        msg = f"failed to locate temp directory: {str(e)}"
        return ResponseData(Status.Failure, message=msg).to_dict()

    # locate file containing included substrates file
    # return error if not found
    try:
        assert os.path.exists(INCLUDED_SUBSTRATES_FILE)
    except Exception as e:
        msg = f"failed to locate included substrates file: {str(e)}"
        return ResponseData(Status.Failure, message=msg).to_dict()

    # parse included substrates files
    # return error if failed
    try:
        included_substrates = parse_substrate_list(INCLUDED_SUBSTRATES_FILE)
        assert included_substrates
    except Exception as e:
        msg = f"Failed to parse included substrates: {str(e)}"
        return ResponseData(Status.Failure, message=msg).to_dict()

    # locate file containing substrate fingerprints file
    # return error if not found
    try:
        assert os.path.exists(FINGERPRINTS_FILE)
    except Exception as e:
        msg = f"failed to locate substrate fingerprints file: {str(e)}"
        return ResponseData(Status.Failure, message=msg).to_dict()

    # parse fingerprints file
    # return error if failed
    try:
        if only_predict_for_smiles_input:
            substrates, fingerprints = [], []  # Only use user-provided SMILES strings.
        else:
            substrates, fingerprints = bitvectors_from_substrate_names(
                substrate_names=included_substrates, path_in_fingerprint_file=FINGERPRINTS_FILE
            )
    except Exception as e:
        msg = f"Failed to parse fingerprints: {str(e)}"
        return ResponseData(Status.Failure, message=msg).to_dict()

    # parse SMILES strings (if provided in frontend)
    # return error if failed
    try:
        if len(smiles_input) > 0:
            smiles_input = smiles_input.strip().split("\n")
            for smiles_string in smiles_input:
                fingerprint = bitvector_from_smiles(
                    smiles=smiles_string, path_in_bitvector_file=FINGERPRINTS_FILE
                )
                fingerprints.append(fingerprint)
                substrates.append(smiles_string)
    except Exception as e:
        msg = f"failed to parse SMILES strings: {str(e)}"
        return ResponseData(Status.Failure, message=msg).to_dict()

    # load model
    # return error if not successful
    try:
        model = joblib.load(os.path.join(MODEL_DIR, "model.parasect"))
        model.set_params(n_jobs=1)
        assert model is not None
    except Exception as e:
        msg = f"failed to load model: {str(e)}"
        return ResponseData(Status.Failure, message=msg).to_dict()

    # run PARASECT
    # return error if not successful
    try:
        results = run_parasect(
            selected_input=selected_input,
            selected_input_type=selected_input_type,
            path_temp_dir=TEMP_DIR,
            model=model,
            substrate_names=substrates,
            substrate_fingerprints=fingerprints,
            use_structure_guided_alignment=use_structure_guided_alignment,
            num_predictions_to_report=num_predictions_to_report,
            first_separator=first_separator,
            second_separator=second_separator,
            third_separator=third_separator,
        )
    except Exception as e:
        msg = f"failed to run PARASECT: {str(e)}"
        return ResponseData(Status.Failure, message=msg).to_dict()

    # return results
    payload = {"results": [r.to_json() for r in results]}
    msg = "submission was successful"

    return ResponseData(Status.Success, message=msg, payload=payload).to_dict()
