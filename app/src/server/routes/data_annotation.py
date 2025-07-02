# -*- coding: utf-8 -*-

"""Routes for making adenylation domain subtrate specificity predictions on raw input."""

import os
import threading
import time
import uuid
from typing import Dict

import joblib
from flask import Blueprint, Response, request, redirect

from parasect.api import run_paras, sort_results

from .app import app
from .common import ResponseData, Status
from .constants import MODEL_DIR, TEMP_DIR

blueprint_annotate_data = Blueprint("annotate_data", __name__)


def run_prediction_protein(job_id: str, data: Dict[str, str]) -> None:
    """Run prediction with PARAS all substrate model on protein input.

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
            selected_input_type = data["selectedInputType"]  # fasta, gbk, or accession
            selected_input = data["selectedInput"]  # fasta src or gbk src

        except Exception as e:
            msg = f"failed to read settings: {str(e)}"
            raise Exception(msg)

        # sanity check selected input type
        # return with message if invalid
        selected_input_type = selected_input_type.strip().lower()
        if selected_input_type not in ["fasta", "gbk", "accession"]:
            msg = f"invalid input type: {selected_input_type}."
            raise Exception(msg)

        # sanity check selected model
        # return with message if invalid

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
            model = joblib.load(os.path.join(MODEL_DIR, "all_substrates_model.paras"))
            model.set_params(n_jobs=1)
            assert model is not None

        except Exception as e:
            msg = f"failed to load model: {str(e)}"
            raise Exception(msg)

        # make predictions
        try:
            results = run_paras(
                selected_input=selected_input,
                selected_input_type=selected_input_type,
                path_temp_dir=TEMP_DIR,
                model=model,
                use_structure_guided_alignment=False,
            )
            for result in results:
                result.sort()
            sorted_results = sort_results(results)

        except Exception as e:
            msg = f"failed to make predictions: {str(e)}"
            raise Exception(msg)

        # clean up, remove loaded model
        del model

        # store results
        new_status = str(Status.Success).lower()
        new_message = "Successfully ran predictions!"
        new_results = [r.to_json() for r in sorted_results]

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


@blueprint_annotate_data.route("/api/annotate_data", methods=["POST"])
def annotate_data() -> Response:
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
    threading.Thread(target=run_prediction_protein, args=(job_id, data)).start()

    # immediately return job_id
    return ResponseData(Status.Success, payload={"jobId": job_id}).to_dict()
