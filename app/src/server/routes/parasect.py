from collections import OrderedDict
from flask import Blueprint, Response, request 
import os
import joblib

from paras.scripts.data_processing.temp import clear_temp
from paras.scripts.feature_extraction.sequence_feature_extraction.seq_to_features import domains_to_features
from paras.scripts.feature_extraction.compound_feature_extraction.fingerprinting import bitvector_from_smiles, bitvectors_from_substrate_names
from paras.scripts.general import get_domains, get_top_n_aa_parasect
from paras.scripts.parsers.parsers import parse_substrate_list

from .common import Status, ResponseData

blueprint_submit_parasect = Blueprint("submit_parasect", __name__)
@blueprint_submit_parasect.route("/api/submit_parasect", methods=["POST"])
def submit_parasect() -> Response:
    """
    Submit settings for prediction with Parasect model.
    
    :return: Response
    """
    data = request.get_json()

    try:
        data = data["data"]
        selected_input = data["src"] # Fasta or Gbk file contents.
        selected_input_type = data["selectedInputType"] # Fasta or Gbk.

        # Options.
        save_active_site_signatures = data["saveActiveSiteSignatures"]
        save_extended_signatures = data["saveExtendedSignatures"]
        save_adenylation_domain_sequences = data["saveAdenylationDomainSequences"]
        num_predictions_to_report = data["numPredictionsToReport"]

        # Advanced options.
        smiles_input = data["smilesSrc"]
        only_predict_for_smiles_input = data["onlyMakePredictionsUploadedSmiles"]
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
        file_name = "input.fasta" if selected_input_type == "fasta" else "input.gbk"
        input_file = os.path.join(temp_dir, file_name)
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
            extraction_method="profile" if use_structure_guided_alignment else "hmm",
            job_name="run",
            separator_1=first_separator,
            separator_2=second_separator,
            separator_3=third_separator,
            verbose=False,
            file_type=selected_input_type.lower(),
            temp_dir=temp_dir
        )
        sequence_ids, sequence_feature_vectors = domains_to_features(a_domains, one_hot=False)

        if not sequence_feature_vectors: 
            raise Exception("No feature vectors.")
        
    except Exception as e:
        clear_temp(temp_dir)
        msg = f"Failed to get domains: {str(e)}"
        return ResponseData(Status.Failure, message=msg).to_dict()
    
    # Locate file containing included substrates.
    try:
        absolute_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        substrate_file = os.path.join(absolute_path, "data/included_substrates.txt")
        assert os.path.exists(substrate_file)

    except Exception as e:
        msg = f"Failed to locate substrate file: {str(e)}"
        return ResponseData(Status.Failure, message=msg).to_dict()
    
    # Parse included substrates.
    try:
        included_substrates = parse_substrate_list(substrate_file)
        assert included_substrates

    except Exception as e:
        msg = f"Failed to parse included substrates: {str(e)}"
        return ResponseData(Status.Failure, message=msg).to_dict()
    
    # Locate file containing substrate fingerprints. 
    try:
        absolute_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        fingerprint_file = os.path.join(absolute_path, "data/fingerprints.txt")
        assert os.path.exists(fingerprint_file)
    
    except Exception as e:
        msg = f"Failed to locate substrate fingerprints: {str(e)}"
        return ResponseData(Status.Failure, message=msg).to_dict()
    
    # Parse fingerprints.
    try:
        if only_predict_for_smiles_input:
            substrates, fingerprints = [], [] # Only use user-provided SMILES strings.
        else:
            substrates, fingerprints = bitvectors_from_substrate_names(included_substrates, fingerprint_file)
    
    except Exception as e:
        msg = f"Failed to parse fingerprints: {str(e)}"
        return ResponseData(Status.Failure, message=msg).to_dict()
    
    # Parse SMILES strings, if provided.
    try:
        if len(smiles_input) > 0:
            smiles_input = smiles_input.strip().split("\n")
            for smiles_string in smiles_input:
                fingerprint = bitvector_from_smiles(smiles_string, fingerprints)
                fingerprints.append(fingerprint)
                substrates.append(smiles_string)
    
    except Exception as e:
        msg = f"Failed to parse SMILES strings: {str(e)}"
        return ResponseData(Status.Failure, message=msg).to_dict()
    
    # Run model and retrieve class predictions.
    try:
        results = OrderedDict()
        if fingerprints and sequence_feature_vectors:

            absolute_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            classifier = joblib.load(os.path.join(absolute_path, "models/model.parasect"))
            classifier.set_params(n_jobs=1)

            batch_size = 1000
            counter = 0
            start = 0
            end = batch_size

            id_to_probabilities = {}

            batch_nr = 1
            while start < len(sequence_feature_vectors):

                labels, feature_vectors = [], []
                for i, sequence_feature_vector in enumerate(sequence_feature_vectors[start:end]):
                    counter += 1
                    for j, fingerprint in enumerate(fingerprints):
                        feature_vector = sequence_feature_vector + fingerprint
                        label = (sequence_ids[start + i], substrates[j])
                        labels.append(label)
                        feature_vectors.append(feature_vector)

                start = counter
                end = min([counter + batch_size, len(sequence_feature_vectors)])

                probabilities = classifier.predict_proba(feature_vectors)
                interaction_labels = classifier.classes_

                if interaction_labels[0] == 1:
                    interaction_index = 0
                elif interaction_labels[1] == 1:
                    interaction_index = 1
                else:
                    raise ValueError("Interaction values must be 0 and 1")

                for i, label in enumerate(labels):
                    seq_id, substrate = label
                    if seq_id not in id_to_probabilities:
                        id_to_probabilities[seq_id] = []

                    id_to_probabilities[seq_id].append((probabilities[i][interaction_index], substrate))

                batch_nr += 1

            for seq_id in sequence_ids:
                results[seq_id] = get_top_n_aa_parasect(seq_id, id_to_probabilities, num_predictions_to_report)

        else:
            classifier = None
            clear_temp()
            msg = "No feature vectors or fingerprints."
            return ResponseData(Status.Failure, message=msg).to_dict()
        
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