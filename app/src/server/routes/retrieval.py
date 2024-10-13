# -*- coding: utf-8 -*-

"""Module for retrieving results from the server."""

from flask import Blueprint

from .app import app
from .common import ResponseData, Status

blueprint_retrieve = Blueprint("retrieve", __name__)


@blueprint_retrieve.route("/api/retrieve/<job_id>", methods=["GET"])
def retrieve(job_id: str):
    """Retrieve results for a job.

    :param job_id: Job ID.
    :type job_id: str
    :return: Response.
    :rtype: Response
    """
    try:
        results = app.config["JOB_RESULTS"][job_id]
        if results["status"] == "pending":
            return ResponseData(Status.Pending, message="job is pending").to_dict()
        elif results["status"] == "failure":
            return ResponseData(Status.Failure, message=results["message"]).to_dict()
        else:
            return ResponseData(Status.Success, payload=results).to_dict()
    except KeyError:
        return ResponseData(Status.Failure, message="job not found").to_dict()
    except Exception as e:
        return ResponseData(Status.Failure, message=str(e)).to_dict()
