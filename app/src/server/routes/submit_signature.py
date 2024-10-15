# -*- coding: utf-8 -*-

"""Routes for making adenylation domain subtrate specificity predictions on signatures."""

import os
import threading
import uuid
from typing import Dict

import joblib
from flask import Blueprint, Response, request

from parasect.api import run_paras_for_signatures
from parasect.core.domain import AdenylationDomain

from .app import app
from .common import ResponseData, Status
from .constants import MODEL_DIR

blueprint_submit_signature = Blueprint("submit_signature", __name__)


def run_prediction_signature(job_id: str, data: Dict[str, str]) -> None:
    """Run prediction with PARAS or PARASECT on signatures.

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
            submissions = data["submissions"]
        except Exception as e:
            msg = f"failed to read settings: {str(e)}"
            raise Exception(msg)

        # sanity checks
        if not isinstance(submissions, list):
            raise Exception("domains must be a list")

        if len(submissions) == 0:
            raise Exception("no domains provided")

        # check that every submission is a dictionary that has at least the keys
        # 'protein_name', 'domain_start', 'domain_end', and 'extended_signature'
        for s in submissions:
            if not isinstance(s, dict):
                raise Exception("each submissions must be a dictionary")
            if not all(
                k in s
                for k in [
                    "protein_name",
                    "domain_start",
                    "domain_end",
                    "signature",
                    "extended_signature",
                ]
            ):
                raise Exception(
                    (
                        "each domain must have keys 'protein_name', 'domain_start', "
                        "'domain_end', 'signature', and 'extended_signature'"
                    )
                )

        # check that every protein_name is a string
        if not all(isinstance(s["protein_name"], str) for s in submissions):
            raise Exception("protein_name must be a string")

        # check that every domain_start is an integer
        if not all(isinstance(s["domain_start"], int) for s in submissions):
            raise Exception("domain_start must be an integer")

        # check that every domain_end is an integer
        if not all(isinstance(s["domain_end"], int) for s in submissions):
            raise Exception("domain_end must be an integer")

        # check that every signature is a string
        if not all(isinstance(s["signature"], str) for s in submissions):
            raise Exception("signature must be a string")

        # check that every extended_signature is a string
        if not all(isinstance(s["extended_signature"], str) for s in submissions):
            raise Exception("extended_signature must be a string")

        # load model
        # return error if not successful
        try:
            model = joblib.load(os.path.join(MODEL_DIR, "model.paras"))
            model.set_params(n_jobs=1)
            assert model is not None
        except Exception as e:
            msg = f"failed to load model: {str(e)}"
            raise Exception(msg)

        # make predictions
        try:
            # create AdenylationDomain objects
            domains = []
            for s in submissions:
                domain = AdenylationDomain(
                    protein_name=s["protein_name"],
                    domain_start=s["domain_start"],
                    domain_end=s["domain_end"],
                )
                domain.signature = s["signature"]
                domain.extended_signature = s["extended_signature"]
                domains.append(domain)

            # sort domains by protein_name and domain_start
            domains = sorted(domains, key=lambda d: (d.protein_name, d.start))

            # asign domain_nr to each domain, count per protein_name
            domain_nr = 0
            protein_name = None
            for domain in domains:
                if domain.protein_name != protein_name:
                    domain_nr = 0
                    protein_name = domain.protein_name
                domain_nr += 1
                domain.set_domain_number(domain_nr)

            # run predictions
            results = run_paras_for_signatures(domains=domains, model=model)

        except Exception as e:
            msg = f"failed to make predictions: {str(e)}"
            raise Exception(msg)

        # clean up, remove loaded model
        del model

        # store results
        app.config["JOB_RESULTS"][job_id] = {
            "status": str(Status.Success).lower(),
            "message": "successfully ran predictions",
            "results": [r.to_json() for r in results],
        }

    except Exception as e:
        app.config["JOB_RESULTS"][job_id] = {
            "status": str(Status.Failure).lower(),
            "message": str(e),
            "results": [],
        }


@blueprint_submit_signature.route("/api/submit_signature", methods=["POST"])
def submit_signature() -> Response:
    """Submit settings for predicting adenylation domain substrate specificities.

    :return: Response.
    :rtype: Response
    """
    data = request.get_json()

    # generate job id
    job_id = str(uuid.uuid4())

    # initialize job  with status as pending
    app.config["JOB_RESULTS"][job_id] = {
        "status": str(Status.Pending).lower(),
        "message": "job is pending",
        "results": [],
    }

    # run prediction in a separate thread
    threading.Thread(target=run_prediction_signature, args=(job_id, data)).start()

    # immediately return job id
    return ResponseData(Status.Success, payload={"jobId": job_id}).to_dict()
