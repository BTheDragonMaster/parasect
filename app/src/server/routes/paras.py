from collections import OrderedDict
from flask import Blueprint, Response, request 
import joblib 
import os

from parasect.api import run_paras

from .common import Status, ResponseData

blueprint_submit_paras = Blueprint("submit_paras", __name__)
@blueprint_submit_paras.route("/api/submit_paras", methods=["POST"])
def submit_paras() -> Response:
    """
    Submit settings for prediction with Paras model.
    
    :return: Response
    """
    data = request.get_json()

    # Read settings. Is everything present?
    try:
        data = data["data"]
        selected_input = data["src"] # Fasta or Gbk file contents.
        selected_input_type = data["selectedInputType"] # Fasta or Gbk.

        # Options.
        save_active_site_signatures = data["saveActiveSiteSignatures"]
        save_extended_signatures = data["saveExtendedSignatures"]
        save_adenylation_domain_sequences = data["saveAdenylationDomainSequences"]
        selected_model = data["selectedSubstrateChoice"] # allSubstrates or commonSubstrates.
        num_predictions_to_report = data["numPredictionsToReport"]

        # Advanced options.
        use_structure_guided_alignment = data["useStructureGuidedProfileAlignment"]
        first_separator = data["firstSeparator"]
        second_separator = data["secondSeparator"]
        third_separator = data["thirdSeparator"]

    except Exception as e:
        msg = f"Failed to read settings: {str(e)}"
        return ResponseData(Status.Failure, message=msg).to_dict()
    
    # Sanity check selected input type.
    selected_input_type = selected_input_type.strip().lower()
    if selected_input_type not in ["fasta", "gbk"]:
        msg = f"Invalid input type: {selected_input_type}."
        return ResponseData(Status.Failure, message=msg).to_dict()
    
    # Sanity check selected model.
    selected_model = selected_model.strip()
    if selected_model not in ["allSubstrates", "commonSubstrates"]:
        msg = f"Invalid model: {selected_model}."
        return ResponseData(Status.Failure, message=msg).to_dict()
    
    # Locate temp directory.
    try:
        absolute_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        temp_dir = os.path.join(absolute_path, "temp")
        assert os.path.exists(temp_dir)

    except Exception as e:
        msg = f"Failed to locate temp directory: {str(e)}"
        return ResponseData(Status.Failure, message=msg).to_dict()
    
    try:
        if selected_model == "commonSubstrates":
            absolute_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            classifier = joblib.load(os.path.join(absolute_path, "models/model.paras"))
            classifier.set_params(n_jobs=1)
            assert classifier is not None

        elif selected_model == "allSubstrates":
            absolute_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            classifier = joblib.load(os.path.join(absolute_path, "models/all_substrates_model.paras"))
            classifier.set_params(n_jobs=1)
            assert classifier is not None
            
        else:
            raise Exception("Invalid model.")

    except Exception as e:
        msg = f"Failed to load model: {str(e)}"
        return ResponseData(Status.Failure, message=msg).to_dict()

    try:
        results = run_paras(
            selected_input=selected_input,
            selected_input_type=selected_input_type,
            temp_dir=temp_dir,
            first_separator=first_separator,
            second_separator=second_separator,
            third_separator=third_separator,
            model=classifier,
            use_structure_guided_alignment=use_structure_guided_alignment,
            num_predictions_to_report=num_predictions_to_report,
            save_active_site_signatures=save_active_site_signatures,
            save_extended_signatures=save_extended_signatures,
            save_adenylation_domain_sequences=save_adenylation_domain_sequences
        )
    except Exception as e:
        msg = f"Failed to run model: {str(e)}"
        return ResponseData(Status.Failure, message=msg).to_dict()
    
    payload = {"results": results}

    msg = "Submission was successful."
    return ResponseData(Status.Success, message=msg, payload=payload).to_dict()