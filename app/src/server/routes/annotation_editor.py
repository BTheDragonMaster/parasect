import logging
from typing import Any
from flask import Blueprint, Response, request, jsonify
from pprint import pprint
import uuid
import os


from parasect.database.query_database import get_substrates_from_smiles, \
    get_substrates_from_name, get_all_substrates, get_protein_names, get_domains_from_synonym
from parasect.database.build_database import Substrate
from pikachu.general import read_smiles

from .database import get_db

blueprint_check_smiles = Blueprint("check_smiles", __name__)
blueprint_check_substrate_name = Blueprint("check_substrate_name", __name__)
blueprint_check_domain_name = Blueprint("check_domain_name", __name__)
blueprint_get_substrates = Blueprint("get_substrates", __name__)
blueprint_submit_annotations = Blueprint("submit_annotations", __name__)
blueprint_check_protein_name = Blueprint("check_protein_name", __name__)


def substrates_from_smiles(smiles: str) -> list[Substrate]:
    session_generator = get_db()
    session = next(session_generator)
    try:
        substrate_hits = get_substrates_from_smiles(session, smiles)
        print(substrate_hits, smiles)
        return substrate_hits

    finally:
        session_generator.close()


def protein_name_in_dataset(protein_name: str) -> bool:
    session_generator = get_db()
    session = next(session_generator)

    try:
        protein_names = get_protein_names(session)
        if protein_name in protein_names:
            return True
        else:
            return False

    finally:
        session_generator.close()


def domain_name_in_dataset(domain_name: str) -> bool:
    session_generator = get_db()
    session = next(session_generator)

    try:
        domains = get_domains_from_synonym(session, domain_name)
        if domains:
            return True
        else:
            return False

    finally:
        session_generator.close()


def substrates_from_name(substrate_name: str) -> list[Substrate]:
    session_generator = get_db()
    session = next(session_generator)
    try:
        substrate_hits = get_substrates_from_name(session, substrate_name)
        return substrate_hits

    finally:
        session_generator.close()


def get_substrate_list() -> list[Substrate]:
    session_generator = get_db()
    session = next(session_generator)
    try:
        substrates = get_all_substrates(session)
        return substrates

    finally:
        session_generator.close()


@blueprint_get_substrates.route("/api/get_substrates", methods=["GET"])
def get_substrates() -> Response:
    substrates = get_substrate_list()
    if not substrates:
        return jsonify({"error": "Missing 'smiles' in request body"}), 400

    json_substrates = [s.to_json() for s in substrates]
    return jsonify(json_substrates)


@blueprint_check_smiles.route("/api/check_smiles", methods=["POST"])
def check_smiles() -> Response:
    """Handle SMILES structure checks and return matching substrates."""
    data = request.get_json()
    if not data or "smiles" not in data:

        return jsonify({"error": "Missing 'smiles' in request body"}), 400
    smiles = data["smiles"]
    try:
        read_smiles(smiles)

    except Exception:

        return jsonify({"invalid": True})

    substrates = substrates_from_smiles(smiles)

    return jsonify([substrate.to_json() for substrate in substrates])


@blueprint_check_substrate_name.route("/api/check_substrate_name", methods=["POST"])
def check_substrate_name() -> Response:
    """Check substrate name against database and return matching substrates."""
    data = request.get_json()
    if not data or "substrate_name" not in data:
        return jsonify({"error": "Missing 'substrate_name' in request body"}), 400
    substrate_name = data["substrate_name"]
    substrates = substrates_from_name(substrate_name)

    return jsonify([substrate.to_json() for substrate in substrates])


@blueprint_check_protein_name.route("/api/check_protein_name", methods=["POST"])
def check_protein_name() -> Response:
    """Check protein name against database """
    data = request.get_json()
    if not data or "protein_name" not in data:
        return jsonify({"error": "Missing 'protein_name' in request body"}), 400
    protein_name = data["protein_name"]
    in_dataset = protein_name_in_dataset(protein_name)

    return jsonify({"protein_in_dataset": in_dataset})


@blueprint_check_domain_name.route("/api/check_domain_name", methods=["POST"])
def check_domain_name() -> Response:
    """Check domain name against database """
    data = request.get_json()
    if not data or "domain_name" not in data:
        return jsonify({"error": "Missing 'domain_name' in request body"}), 400
    domain_name = data["domain_name"]
    in_dataset = domain_name_in_dataset(domain_name)

    return jsonify({"domain_in_dataset": in_dataset})


def cleanup_annotations(annotations: dict[str, Any]):
    protein_to_entries = {}

    for protein in annotations.keys():

        domain_to_annotations = annotations[protein]["domains"]

        synonym = annotations[protein]["synonym"]
        domain_to_entries = {}
        for domain, domain_annotations in domain_to_annotations.items():
            domain_name = domain_annotations["name"]
            for annotation in domain_annotations["substrates"]:
                if 'substrateName' in annotation and 'substrateSmiles' in annotation and\
                        annotation['substrateName'] and annotation['substrateSmiles'] and \
                        annotation['annotationType'] and annotation['sequence'] and \
                        annotation['signature'] and annotation['extendedSignature'] and \
                        annotation['annotationType'] != 'no_update':
                    if domain_name not in domain_to_entries:
                        domain_to_entries[domain_name] = {"sequence": annotation['sequence'],
                                                          "signature": annotation['signature'],
                                                          "extended_signature": annotation['extendedSignature'],
                                                          "substrates": []}
                    domain_to_entries[domain_name]["substrates"].append({"name": annotation['substrateName'],
                                                                         "smiles": annotation['substrateSmiles'],
                                                                         "annotation_type": annotation['annotationType']})

        if domain_to_entries:
            protein_to_entries[synonym] = domain_to_entries

    return protein_to_entries


def write_annotations(file_name, annotations):
    with open(file_name, 'w') as out:
        out.write("domain_id\tsequence\tsignature\textended_signature\tspecificity\tsmiles\tannotation_type\n")
        for protein, domain_annotations in annotations.items():
            for domain_id, domain_info in domain_annotations.items():
                substrates = domain_info["substrates"]

                specificities = []
                smiles = []

                for substrate_info in substrates:
                    annotation_type = substrate_info["annotation_type"]
                    specificities.append(substrate_info["name"])

                    smiles.append(substrate_info["smiles"])

                specificities = '|'.join(specificities)
                smiles = '|'.join(smiles)
                out.write(f"{domain_id}\t{domain_info['sequence']}\t{domain_info['signature']}\t{domain_info['extended_signature']}\t{specificities}\t{smiles}\t{annotation_type}\n")


@blueprint_submit_annotations.route("/api/submit_annotations", methods=["POST"])
def submit_annotations():
    try:
        data = request.get_json()
        annotations = data.get("annotations", {})
        protein_to_entries = cleanup_annotations(annotations)
        pprint(protein_to_entries)
        print("hello!")
        current_dir = os.path.dirname(os.path.abspath(__file__))
        user_submission = os.path.join(current_dir, '..', 'user_submissions', f"{str(uuid.uuid4())}.txt")
        write_annotations(user_submission, protein_to_entries)



        # Future: validate or transform annotations here
        # Example: check structure, check duplicates, etc.

        # Create a GitHub PR with annotations
        pr_url = create_github_pr(annotations)

        return jsonify({"pr_url": pr_url}), 200

    except Exception as e:
        logging.exception("Error in /api/submit_annotations")
        return jsonify({"error": str(e)}), 500


def create_github_pr(annotations):
    """
    TODO: Replace this stub with logic to:
    - write a JSON or YAML file from the annotations
    - commit it to a new branch
    - open a GitHub pull request using PyGithub or Git CLI
    """
    # Simulate PR creation
    print("Annotations received for PR:", annotations)

    # Replace this with actual PR URL
    return "https://github.com/your-org/your-repo/pull/123"


