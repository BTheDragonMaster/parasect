# -*- coding: utf-8 -*-

"""Routes for making adenylation domain subtrate specificity predictions with PARAS."""

import os

import joblib
from flask import Blueprint, Response, request

from parasect.api import run_paras

from .common import ResponseData, Status
from .constants import MODEL_DIR, TEMP_DIR

blueprint_submit_paras = Blueprint("submit_paras", __name__)


@blueprint_submit_paras.route("/api/submit_paras", methods=["POST"])
def submit_paras() -> Response:
    """Submit settings for prediction with Paras model.

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
        selected_model = data["selectedSubstrateChoice"]  # allSubstrates or commonSubstrates.
        num_predictions_to_report = data["numPredictionsToReport"]
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

    # sanity check selected model
    # return with message if invalid
    selected_model = selected_model.strip()
    if selected_model not in ["allSubstrates", "commonSubstrates"]:
        msg = f"invalid model: {selected_model}"
        return ResponseData(Status.Failure, message=msg).to_dict()

    # locate temp directory
    # return error if not found.
    try:
        assert os.path.exists(TEMP_DIR)
    except Exception as e:
        msg = f"Failed to locate temp directory: {str(e)}"
        return ResponseData(Status.Failure, message=msg).to_dict()

    # load model
    # return error if not successful
    try:
        if selected_model == "commonSubstrates":
            model = joblib.load(os.path.join(MODEL_DIR, "model.paras"))
            model.set_params(n_jobs=1)
            assert model is not None
        elif selected_model == "allSubstrates":
            model = joblib.load(os.path.join(MODEL_DIR, "all_substrates_model.paras"))
            model.set_params(n_jobs=1)
            assert model is not None
        else:
            raise Exception("not successful loading model")
    except Exception as e:
        msg = f"Failed to load model: {str(e)}"
        return ResponseData(Status.Failure, message=msg).to_dict()

    # run PARAS
    # return error if not successful
    try:
        results = run_paras(
            selected_input=selected_input,
            selected_input_type=selected_input_type,
            path_temp_dir=TEMP_DIR,
            first_separator=first_separator,
            second_separator=second_separator,
            third_separator=third_separator,
            model=model,
            use_structure_guided_alignment=use_structure_guided_alignment,
            num_predictions_to_report=num_predictions_to_report,
        )
    except Exception as e:
        msg = f"failed to run PARAS: {str(e)}"
        return ResponseData(Status.Failure, message=msg).to_dict()

    # return results
    payload = {"results": [r.to_json() for r in results]}
    msg = "submission was successful"

    return ResponseData(Status.Success, message=msg, payload=payload).to_dict()
