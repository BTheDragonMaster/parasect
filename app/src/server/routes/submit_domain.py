# -*- coding: utf-8 -*-

"""Route for predicting adenylation domain substrate specificity from a full domain sequence.

`/api/submit_quick` (see submit.py) expects the caller to have already
extracted the 34-residue extended signature itself. This route accepts
the full adenylation domain amino acid sequence and mines the signature
server-side, via the same MUSCLE profile-alignment pipeline already used
by the "structure guided alignment" option on `/api/submit_raw`
(`AdenylationDomain.set_domain_signatures_profile`).
"""

import threading
import time
import uuid
from typing import Dict

from flask import Blueprint, Response, jsonify, request

from parasect.api import run_paras_for_signatures
from parasect.core.domain import AdenylationDomain
from parasect.core.hit import DomainType

from .common import ResponseData, Status
from .constants import cleanup_job_temp_dir, job_temp_dir
from .job_store import set_job, update_job
from .submit import loader

blueprint_submit_domain = Blueprint("submit_domain", __name__)

# only whole-domain types make sense for a standalone domain sequence:
# AMP_BINDING_C and N_TERMINAL are HMM hit classifications, not domain types
_DOMAIN_TYPES = {
    "AMP_BINDING": DomainType.AMP_BINDING,
    "A_OX": DomainType.A_OX,
}


def run_prediction_domain(job_id: str, data: Dict[str, str]) -> None:
    """Mine extended signatures from full domain sequences, then predict.

    :param job_id: Job ID.
    :param data: Data.
    """
    path_temp_dir = job_temp_dir(job_id)

    try:
        try:
            data = data["data"]
            submissions = data["submissions"]
        except Exception as e:
            raise Exception(f"failed to read settings: {str(e)}")

        if not isinstance(submissions, list):
            raise Exception("submissions must be a list")

        if len(submissions) == 0:
            raise Exception("no domains provided")

        for s in submissions:
            if not isinstance(s, dict):
                raise Exception("each submission must be a dictionary")
            if not all(k in s for k in ["protein_name", "domain_sequence"]):
                raise Exception("each submission must have keys 'protein_name' and 'domain_sequence'")
            if not isinstance(s["protein_name"], str):
                raise Exception("protein_name must be a string")
            if not isinstance(s["domain_sequence"], str) or not s["domain_sequence"].strip():
                raise Exception("domain_sequence must be a non-empty string")
            domain_type = s.get("domain_type", "AMP_BINDING")
            if domain_type not in _DOMAIN_TYPES:
                raise Exception(
                    f"invalid domain_type '{domain_type}', must be one of {sorted(_DOMAIN_TYPES)}"
                )

        try:
            model = loader.get("parasAllSubstrates")
        except Exception as e:
            raise Exception(f"failed to load model: {str(e)}")

        try:
            domains = []
            for s in submissions:
                sequence = s["domain_sequence"].strip()
                domain = AdenylationDomain(
                    protein_name=s["protein_name"],
                    domain_type=_DOMAIN_TYPES[s.get("domain_type", "AMP_BINDING")],
                    domain_start=0,
                    domain_end=len(sequence),
                )
                domain.set_sequence(sequence)
                domain.set_protein_sequence(sequence)

                try:
                    domain.set_domain_signatures_profile(path_temp_dir)
                except Exception as e:
                    raise Exception(
                        f"failed to mine signature for domain '{s['protein_name']}': {str(e)}"
                    )

                if not domain.extended_signature:
                    raise Exception(
                        f"could not extract a 34-residue signature for domain '{s['protein_name']}'; "
                        "check that the sequence is a valid adenylation domain"
                    )

                domains.append(domain)

            # sort domains by protein_name and domain_start, then assign
            # per-protein domain numbers, same convention as /api/submit_quick
            domains = sorted(domains, key=lambda d: (d.protein_name, d.start))
            domain_nr = 0
            protein_name = None
            for domain in domains:
                if domain.protein_name != protein_name:
                    domain_nr = 0
                    protein_name = domain.protein_name
                domain_nr += 1
                domain.set_domain_number(domain_nr)

            results = run_paras_for_signatures(domains=domains, model=model)

        except Exception as e:
            raise Exception(f"failed to make predictions: {str(e)}")

        del model

        update_job(
            job_id,
            status=str(Status.Success).lower(),
            message="Successfully ran predictions!",
            results=[r.to_json() for r in results],
        )

    except Exception as e:
        update_job(job_id, status=str(Status.Failure).lower(), message=str(e), results=[])

    finally:
        cleanup_job_temp_dir(job_id)


@blueprint_submit_domain.route("/api/submit_domain", methods=["POST"])
def submit_domain() -> Response:
    """Submit full adenylation domain sequences for signature mining + prediction.

    Expects a JSON body shaped like:
    ``{"data": {"submissions": [{"protein_name": str, "domain_sequence": str,
    "domain_type": "AMP_BINDING" | "A_OX" (optional, default "AMP_BINDING")}, ...]}}``

    Unlike `/api/submit_quick`, callers do not need to extract the 34-residue
    extended signature themselves -- it's mined server-side via MUSCLE
    profile alignment. Async, same contract as `/api/submit_raw`: returns a
    `jobId` immediately, poll `/api/retrieve/<job_id>` for the result.

    :return: Response.
    """
    data = request.get_json()

    job_id = str(uuid.uuid4())
    current_time = int(time.time())

    set_job(job_id, {
        "status": str(Status.Pending).lower(),
        "message": "Job is pending!",
        "results": [],
        "timestamp": current_time,
    })

    threading.Thread(target=run_prediction_domain, args=(job_id, data)).start()

    return jsonify(ResponseData(Status.Success, payload={"jobId": job_id}).to_dict())
