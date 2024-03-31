from collections import OrderedDict
from flask import Blueprint, Response, request 
import joblib 
import os

from paras.scripts.data_processing.temp import clear_temp
from paras.scripts.feature_extraction.sequence_feature_extraction.seq_to_features import domains_to_features
from paras.scripts.general import get_domains, get_top_n_aa_paras

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

    except Exception as e:
        msg = f"Failed to read settings: {str(e)}"
        return ResponseData(Status.Failure, message=msg).to_dict()
    
    # Sanity check selected input type.
    selected_input_type = selected_input_type.strip().lower()
    if selected_input_type not in ["fasta", "genbank"]:
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

    # Try to write selected_input to file in temp folder.
    try:
        input_file = os.path.join(temp_dir, "input.fasta")
        with open(input_file, "w") as f:
            f.write(selected_input)

    except Exception as e:
        clear_temp(temp_dir)
        msg = f"Failed to write input to file: {str(e)}"
        return ResponseData(Status.Failure, message=msg).to_dict()
    
    # Get domains.
    try:
        a_domains = get_domains(
            input_file=input_file, 
            extraction_method="hmm", 
            job_name="run", 
            separator_1=first_separator, 
            separator_2=second_separator, 
            separator_3=third_separator, 
            verbose=False,
            file_type=selected_input_type.lower(), 
            temp_dir=temp_dir
        )
        ids, feature_vectors = domains_to_features(a_domains, one_hot=False)
        if not feature_vectors: 
            raise Exception("No feature vectors.")
        
    except Exception as e:
        clear_temp(temp_dir)
        msg = f"Failed to get domains: {str(e)}"
        return ResponseData(Status.Failure, message=msg).to_dict()
    
    # Run model and retrieve class predictions.
    try:
        results = OrderedDict()
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

        probabilities = classifier.predict_proba(feature_vectors)
        amino_acid_classes = classifier.classes_

        for i, seq_id in enumerate(ids):
            probability_list = probabilities[i]
            probs_and_aa = get_top_n_aa_paras(amino_acid_classes, probability_list, num_predictions_to_report)
            results[seq_id] = probs_and_aa

        classifier = None
        clear_temp(temp_dir)

    except Exception as e:
        classifier = None
        clear_temp(temp_dir)
        msg = f"Failed to run model: {str(e)}"
        return ResponseData(Status.Failure, message=msg).to_dict()
    
    # Parse results
    try:
        domain_results = {}
        for domain in a_domains:
            domain_results[domain.domain_id] = {}
            if save_adenylation_domain_sequences:
                domain_results[domain.domain_id]['sequence'] = domain.sequence
            if save_active_site_signatures:
                domain_results[domain.domain_id]['signature'] = domain.signature
            if save_extended_signatures:
                domain_results[domain.domain_id]['extended_signature'] = domain.extended_signature

        for domain_id in results:
            preds = results[domain_id]
            domain_results[domain_id]['predictions'] = preds

    except Exception as e:
        msg = f"Failed to parse results: {str(e)}"
        return ResponseData(Status.Failure, message=msg).to_dict()
    
    # Create payload.
    payload = {"results": [
        {
            "domain_id": domain_id,
            "data": domain_results[domain_id]
        } for domain_id in domain_results
    ]}

    msg = "Submission was successful."
    return ResponseData(Status.Success, message=msg, payload=payload).to_dict()