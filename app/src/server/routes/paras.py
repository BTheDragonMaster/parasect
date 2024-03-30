from flask import Blueprint, Response, request 
import time

from .common import Status, ResponseData

blueprint_submit_paras = Blueprint("submit_paras", __name__)
@blueprint_submit_paras.route("/api/submit_paras", methods=["POST"])
def submit_paras() -> Response:
    """
    Submit settings for prediction with Paras model.
    
    :return: Response
    """
    data = request.get_json()

    try:
        data = data["data"]
        selected_input = data["src"] # Fasta or Genbank file contents.
        selected_input_type = data["selectedInputType"] # Fasta or Genbank.

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

    except:
        msg = "Failed to read settings."
        return ResponseData(Status.Error, message=msg).to_dict()

    msg = "Submission was successful."
    return ResponseData(Status.Success, message=msg).to_dict()