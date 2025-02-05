# -*- coding: utf-8 -*-

"""Routes for making adenylation domain subtrate specificity predictions on raw input."""

import os
import threading
import time
import uuid
from typing import Dict

import joblib
from flask import Blueprint, Response, request, redirect

from parasect.api import run_paras, run_parasect, run_paras_for_signatures
from parasect.core.domain import AdenylationDomain
from pikachu.general import read_smiles

from .app import app
from .common import ResponseData, Status
from .constants import MODEL_DIR, TEMP_DIR

blueprint_submit_raw = Blueprint("submit_raw", __name__)
blueprint_submit_quick = Blueprint("submit_quick", __name__)


########################################################################################################################
########################################################################################################################
#
# Submit settings with JSON parameters
#
########################################################################################################################
########################################################################################################################


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
            uploaded_substrates_file_has_header = data["uploadedSubstratesFileContentHasHeader"]
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
                        for line_index, line in enumerate(lines):
                            
                            # skip header line if present
                            if uploaded_substrates_file_has_header and line_index == 0:
                                continue

                            smiles_name, smiles = line.strip().split("\t")

                            # check if SMILES string is valid
                            try:
                                _ = read_smiles(smiles)
                            except Exception as e:
                                msg = f"failed to parse SMILES string from custom substrates on line {line_index + 1} with name '{smiles_name}'. Error: {str(e)}"
                                raise Exception(msg)

                            custom_substrate_names.append(smiles_name)
                            custom_substrate_smiles.append(smiles)
                except Exception as e:
                    msg = f"failed to parse SMILES strings: {str(e)}"
                    raise Exception(msg)
                
                # if no custom substrates were provided, and also the option to use only uploaded substrates is selected
                # return with error message that at least one custom substrate must be provided in that case
                if len(custom_substrate_names) == 0 and use_only_uploaded_substrates:
                    raise Exception("at least one custom substrate must be provided if 'Use only uploaded substrates' is selected.")

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


########################################################################################################################
########################################################################################################################
#
# Submit settings with signature parameters directly in URL
#
########################################################################################################################
########################################################################################################################


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
                    "extended_signature",
                ]
            ):
                raise Exception(
                    (
                        "each domain must have keys 'protein_name', 'domain_start', "
                        "'domain_end', and 'extended_signature'"
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


@blueprint_submit_quick.route("/api/submit_quick", methods=["GET"])
def submit_quick() -> Response:
    """Submit settings with signature parameters directly in URL and redirect to results page.
    
    Example URL for localhost, if backend is running on port 5000:
    http://localhost:5000/api/submit_quick?signature1=LDQIFDVFVSEMSLIVGGEVNAYGPTETTVEATA&name1=NameA&start1=10&end1=50

    Example URL for production:
    https://paras.bioinformatics.nl/api/submit_quick?signature1=LDQIFDVFVSEMSLIVGGEVNAYGPTETTVEATA&name1=NameA&start1=10&end1=50
    """        
    job_id = str(uuid.uuid4())

    try:
        signature_keys = [key for key in request.args.keys() if key.startswith("signature")]

        # parse submissions from URL parameters
        data = {"data": {"submissions": []}}
        for signature_key in signature_keys:
            # extract accession from the signature key
            accession = signature_key.replace("signature", "")

            # retrieve actual signature value
            signature = request.args.get(signature_key, None)

            if signature is None:
                return f"no signature provided for {accession}", 400

            # retrieve optional attributes based on accession
            name = request.args.get(f"name{accession}", f"Signature {accession}")
            start = request.args.get(f"start{accession}", 0, type=int)
            end = request.args.get(f"end{accession}", 0, type=int)
            submission = {
                "protein_name": name,
                "domain_start": start,
                "domain_end": end,
                "extended_signature": signature,
            }
            data["data"]["submissions"].append(submission)

        # sort submissions first by domain start, then by protein name
        data["data"]["submissions"] = sorted(
            data["data"]["submissions"], 
            key=lambda x: (x["domain_start"], x["protein_name"])
        )

        if len(data["data"]["submissions"]) == 0:
            app.config["JOB_RESULTS"][job_id] = {
                "status": str(Status.Failure).lower(),
                "message": "No valid signatures provided.",
                "results": [],
                "timestamp": int(time.time()),
            }

        else:
            # Initialize job with pending status
            app.config["JOB_RESULTS"][job_id] = {
                "status": str(Status.Pending).lower(),
                "message": "Job is pending!",
                "results": [],
                "timestamp": int(time.time()),
            }

            threading.Thread(target=run_prediction_signature, args=(job_id, data)).start()
    
    except Exception as e:
        app.config["JOB_RESULTS"][job_id] = {
            "status": str(Status.Failure).lower(),
            "message": str(e),
            "results": [],
            "timestamp": int(time.time()),
        }

    if app.config["ENV"] == "production":
        results_url = f"https://paras.bioinformatics.nl/results/{job_id}"
    else:
        # frontend is running on localhost port 3000
        results_url = f"http://localhost:3000/results/{job_id}"

    return redirect(results_url)

        
