import os
from collections import OrderedDict

import joblib
from flask import Blueprint, Response, request

from parasect.api import run_parasect
from parasect.core.featurisation import (
    bitvector_from_smiles,
    bitvectors_from_substrate_names,
)
from parasect.core.helpers import clear_temp_dir
from parasect.core.parsing import parse_substrate_list

from .common import ResponseData, Status

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
        selected_input = data["src"]  # Fasta or Gbk file contents.
        selected_input_type = data["selectedInputType"]  # Fasta or Gbk.

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
        clear_temp_dir(temp_dir, keep=[".gitkeep"])
        msg = f"Failed to write input to file: {str(e)}"
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
            substrates, fingerprints = [], []  # Only use user-provided SMILES strings.
        else:
            substrates, fingerprints = bitvectors_from_substrate_names(
                included_substrates, fingerprint_file
            )

    except Exception as e:
        msg = f"Failed to parse fingerprints: {str(e)}"
        return ResponseData(Status.Failure, message=msg).to_dict()

    # Parse SMILES strings, if provided.
    try:
        if len(smiles_input) > 0:
            smiles_input = smiles_input.strip().split("\n")
            for smiles_string in smiles_input:
                bitvector_file_path = os.path.join(absolute_path, "data/fingerprints.txt")
                fingerprint = bitvector_from_smiles(smiles_string, bitvector_file_path)
                fingerprints.append(fingerprint)
                substrates.append(smiles_string)

    except Exception as e:
        msg = f"Failed to parse SMILES strings: {str(e)}"
        return ResponseData(Status.Failure, message=msg).to_dict()

    try:
        absolute_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        classifier = joblib.load(os.path.join(absolute_path, "models/model.parasect"))
        classifier.set_params(n_jobs=1)
        assert classifier is not None

    except Exception as e:
        msg = f"Failed to load model: {str(e)}"
        return ResponseData(Status.Failure, message=msg).to_dict()

    try:
        results = run_parasect(
            model=classifier,
            input_file=input_file,
            selected_input_type=selected_input_type,
            first_separator=first_separator,
            second_separator=second_separator,
            third_separator=third_separator,
            use_structure_guided_alignment=use_structure_guided_alignment,
            fingerprints=fingerprints,
            temp_dir=temp_dir,
            substrates=substrates,
            num_predictions_to_report=num_predictions_to_report,
            save_active_site_signatures=save_active_site_signatures,
            save_extended_signatures=save_extended_signatures,
            save_adenylation_domain_sequences=save_adenylation_domain_sequences,
        )
    except Exception as e:
        msg = f"Failed to run model: {str(e)}"
        return ResponseData(Status.Failure, message=msg).to_dict()

    payload = {"results": results}

    msg = "Submission was successful."
    return ResponseData(Status.Success, message=msg, payload=payload).to_dict()
