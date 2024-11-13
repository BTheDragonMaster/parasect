# -*- coding: utf-8 -*-

"""Routes for making adenylation domain subtrate specificity predictions on raw input."""

import os
import threading
import time
import uuid
from typing import Dict

import joblib
from flask import Blueprint, Response, request

from parasect.api import run_paras, run_parasect

from .app import app
from .common import ResponseData, Status
from .constants import MODEL_DIR, TEMP_DIR

blueprint_submit_raw = Blueprint("submit_raw", __name__)
blueprint_submit_quick = Blueprint("submit_quick", __name__)


def run_prediction_raw(job_id: str, data: Dict[str, str]) -> None:
    """Run prediction with PARAS or PARASECT on raw data.

    :param job_id: Job ID.
    :type job_id: str
    :param data: Data.
    :type data: Dict[str, str]
    """
    try:
        # read settings
        # is everything present? return with message if not
        try:
            data = data["data"]
            selected_input_type = data["selectedInputType"]  # fasta or gbk
            selected_input = data["selectedInput"]  # fasta src or gbk src
            selected_model = data[
                "selectedModel"
            ]  # parasAllSubstrates, parasCommonSubstrates, or parasect
            use_structure_guided_alignment = data["useStructureGuidedAlignment"]
            smiles_file_content = data["smilesFileContent"]
            use_only_uploaded_substrates = data["useOnlyUploadedSubstrates"]
        except Exception as e:
            msg = f"failed to read settings: {str(e)}"
            raise Exception(msg)

        # sanity check selected input type
        # return with message if invalid
        selected_input_type = selected_input_type.strip().lower()
        if selected_input_type not in ["fasta", "gbk"]:
            msg = f"invalid input type: {selected_input_type}."
            raise Exception(msg)

        # sanity check selected model
        # return with message if invalid
        selected_model = selected_model.strip()
        if selected_model not in ["parasAllSubstrates", "parasCommonSubstrates", "parasect"]:
            msg = f"invalid model: {selected_model}"
            raise Exception(msg)

        # locate temp directory
        # return error if not found.
        try:
            assert os.path.exists(TEMP_DIR)
        except Exception as e:
            msg = f"failed to locate temp directory: {str(e)}"
            raise Exception(msg)

        # load model
        # return error if not successful
        try:
            if selected_model == "parasAllSubstrates":
                model = joblib.load(os.path.join(MODEL_DIR, "all_substrates_model.paras"))
                model.set_params(n_jobs=1)
                assert model is not None
            elif selected_model == "parasCommonSubstrates":
                model = joblib.load(os.path.join(MODEL_DIR, "model.paras"))
                model.set_params(n_jobs=1)
                assert model is not None
            elif selected_model == "parasect":
                model = joblib.load(os.path.join(MODEL_DIR, "model.parasect"))
                model.set_params(n_jobs=1)
                assert model is not None
            else:
                raise Exception("not successful loading model")
        except Exception as e:
            msg = f"failed to load model: {str(e)}"
            raise Exception(msg)

        # make predictions
        try:
            if selected_model in ["parasAllSubstrates", "parasCommonSubstrates"]:

                # run PARAS
                results = run_paras(
                    selected_input=selected_input,
                    selected_input_type=selected_input_type,
                    path_temp_dir=TEMP_DIR,
                    model=model,
                    use_structure_guided_alignment=use_structure_guided_alignment,
                )

            elif selected_model == "parasect":
                # parse user-provided SMILES strings (if provided in frontend)
                # return error if failed
                custom_substrate_names = []
                custom_substrate_smiles = []

                try:
                    if len(smiles_file_content) > 0:
                        lines = smiles_file_content.strip().split("\n")
                        for line in lines:
                            smiles_name, smiles = line.strip().split("\t")
                            custom_substrate_names.append(smiles_name)
                            custom_substrate_smiles.append(smiles)
                except Exception as e:
                    msg = f"failed to parse SMILES strings: {str(e)}"
                    raise Exception(msg)

                # run PARASECT
                results = run_parasect(
                    selected_input=selected_input,
                    selected_input_type=selected_input_type,
                    path_temp_dir=TEMP_DIR,
                    model=model,
                    custom_substrate_names=custom_substrate_names,
                    custom_substrate_smiles=custom_substrate_smiles,
                    only_custom=use_only_uploaded_substrates,
                    use_structure_guided_alignment=use_structure_guided_alignment,
                )

            else:
                raise Exception(f"invalid model {selected_model}")

        except Exception as e:
            msg = f"failed to make predictions: {str(e)}"
            raise Exception(msg)

        # clean up, remove loaded model
        del model

        # store results
        new_status = str(Status.Success).lower()
        new_message = "Successfully ran predictions!"
        new_results = [r.to_json() for r in results]

        app.config["JOB_RESULTS"][job_id]["status"] = new_status
        app.config["JOB_RESULTS"][job_id]["message"] = new_message
        app.config["JOB_RESULTS"][job_id]["results"] = new_results

    except Exception as e:
        # store results
        new_status = str(Status.Failure).lower()
        new_message = str(e)
        new_results = []

        app.config["JOB_RESULTS"][job_id]["status"] = new_status
        app.config["JOB_RESULTS"][job_id]["message"] = new_message
        app.config["JOB_RESULTS"][job_id]["results"] = new_results


@blueprint_submit_raw.route("/api/submit_raw", methods=["POST"])
def submit_raw() -> Response:
    """Submit settings for predicting adenylation domain substrate specificities.

    :return: Response.
    :rtype: Response
    """
    data = request.get_json()

    # generate a uuid
    job_id = str(uuid.uuid4())

    # get current time in seconds since Unix epoch
    current_time = int(time.time())

    # initialize job with status as pending
    app.config["JOB_RESULTS"][job_id] = {
        "status": str(Status.Pending).lower(),
        "message": "Job is pending!",
        "results": [],
        "timestamp": current_time,
    }

    # run prediction in a separate thread
    threading.Thread(target=run_prediction_raw, args=(job_id, data)).start()

    # immediately return job_id
    return ResponseData(Status.Success, payload={"jobId": job_id}).to_dict()


@blueprint_submit_quick.route("/api/submit_quick", methods=["POST"])
def submit_quick() -> Response:
    """Submit settings for predicting adenylation domain substrate specificities.
    
    :return: Response.
    :rtype: Response
    """
    # TODO: implement
    pass
