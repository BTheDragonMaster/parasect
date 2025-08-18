import logging
from typing import Any
from flask import Blueprint, Response, request, jsonify
import uuid
import os


from parasect.database.query_database import get_substrates_from_smiles, \
    get_substrates_from_name, get_all_substrates, get_protein_names, get_domains_from_synonym
from parasect.database.build_database import Substrate
from parasect.core.github import submit_github_issues
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


def cleanup_annotations(annotations: dict[str, dict[str, Any]]) -> dict[str, dict[str, Any]]:
    protein_to_entries: dict[str, dict[str, Any]] = {}

    for protein in annotations.keys():

        domain_to_annotations = annotations[protein]["domains"]

        synonym = annotations[protein]["synonym"]
        domain_to_entries: dict[str, dict[str, Any]] = {}
        for domain, domain_annotations in domain_to_annotations.items():
            domain_name = domain_annotations["name"]
            annotation_type = domain_annotations["annotationType"]
            if annotation_type == 'no_update':
                continue
            elif not annotation_type:
                continue
            else:
                for annotation in domain_annotations["substrates"]:
                    if 'substrateName' in annotation and 'substrateSmiles' in annotation and\
                            annotation['substrateName'] and annotation['substrateSmiles'] and annotation['sequence'] and \
                            annotation['signature'] and annotation['extendedSignature']:
                        if domain_name not in domain_to_entries:
                            domain_to_entries[domain_name] = {"sequence": annotation['sequence'],
                                                              "signature": annotation['signature'],
                                                              "extended_signature": annotation['extendedSignature'],
                                                              "substrates": [],
                                                              "annotation_type": annotation_type}
                        domain_to_entries[domain_name]["substrates"].append({"name": annotation['substrateName'],
                                                                             "smiles": annotation['substrateSmiles'],})

        if domain_to_entries:
            protein_to_entries[synonym] = {"domains": domain_to_entries,
                                           "sequence": annotations[protein]["sequence"]}

    return protein_to_entries


def protein_has_annotation_type(domain_annotations, annotation_name):
    for domain_info in domain_annotations.values():
        if domain_info["annotation_type"] == annotation_name:
            return True

    return False


def data_has_annotation_type(annotations, annotation_name):
    for protein_annotations in annotations.values():
        if protein_has_annotation_type(protein_annotations["domains"], annotation_name):
            return True
    return False


def write_annotations(submission_dir, annotations):
    new_entries = os.path.join(submission_dir, "new")
    corrections = os.path.join(submission_dir, "corrections")
    duplicates = os.path.join(submission_dir, "duplicates")

    all_annotation_types = ["new_entry", "correction", "duplicate_entry"]
    directories = [new_entries, corrections, duplicates]

    annotation_types = []
    submission_dirs = []

    for i, annotation_type in enumerate(all_annotation_types):

        if data_has_annotation_type(annotations, annotation_type):
            submission_dirs.append(directories[i])
            annotation_types.append(annotation_type)

    substrates_out = os.path.join(submission_dir, "smiles.tsv")
    new_substrates = []

    for i, directory in enumerate(submission_dirs):
        annotation_type = annotation_types[i]
        os.mkdir(directory)
        domain_out = os.path.join(directory, "domains.fasta")
        protein_out = os.path.join(directory, "proteins.fasta")
        parasect_data_out = os.path.join(directory, "parasect_data.txt")
        signatures_out = os.path.join(directory, "signatures.fasta")
        extended_signatures_out = os.path.join(directory, "extended_signatures.fasta")
        with open(domain_out, 'w') as d_out, open(protein_out, 'w') as p_out, \
                open(parasect_data_out, 'w') as data_out, open(signatures_out, 'w') as sig_out, \
                open(extended_signatures_out, 'w') as ext_out:
            data_out.write("domain_id\tsequence\tspecificity\n")
            for protein, protein_annotations in annotations.items():
                if not protein_has_annotation_type(protein_annotations["domains"], annotation_type):
                    continue

                p_out.write(f">{protein}\n{protein_annotations['sequence']}\n")
                domain_annotations = protein_annotations["domains"]
                for domain_id, domain_info in domain_annotations.items():
                    if domain_info["annotation_type"] != annotation_type:
                        continue

                    substrates = domain_info["substrates"]
                    specificities = []

                    for substrate_info in substrates:
                        specificities.append(substrate_info["name"])
                        if not substrates_from_name(substrate_info["name"]):
                            new_substrates.append((substrate_info["name"], substrate_info["smiles"]))

                    d_out.write(f">{domain_id}\n{domain_info['sequence']}\n")
                    sig_out.write(f">{domain_id}\n{domain_info['signature']}\n")
                    ext_out.write(f">{domain_id}\n{domain_info['extended_signature']}\n")

                    specificities = '|'.join(specificities)
                    data_out.write(f"{domain_id}\t{domain_info['sequence']}\t{specificities}\n")

    if new_substrates:
        new_substrates = list(set(new_substrates))
        with open(substrates_out, 'w') as s_out:
            s_out.write("substrate\tsmiles\n")
            for new_substrate in new_substrates:
                s_out.write(f"{new_substrate[0]}\t{new_substrate[1]}\n")


@blueprint_submit_annotations.route("/api/submit_annotations", methods=["POST"])
def submit_annotations():
    try:
        data = request.get_json()
        annotations = data.get("annotations", {})
        protein_to_entries = cleanup_annotations(annotations)

        current_dir = os.path.dirname(os.path.abspath(__file__))
        job_id = str(uuid.uuid4())
        user_submission = os.path.join(current_dir, '..', 'user_submissions', job_id)
        os.mkdir(user_submission)
        write_annotations(user_submission, protein_to_entries)

        # Create a GitHub PR with annotations
        pr_url = create_github_issue(protein_to_entries)

        return jsonify({"pr_url": pr_url}), 200

    except Exception as e:
        logging.exception("Error in /api/submit_annotations")
        return jsonify({"error": str(e)}), 500


def create_github_issue(annotations: dict[str, dict[str, Any]]) -> str:
    """
    Create one GitHub issue per protein

    :param annotations: dictionary of protein entries
    :type annotations: dict[str, dict[str, Any]]
    :return: Link to github issues page
    :rtype: str
    """
    print("Annotations received for GitHub issue:", annotations)
    session_generator = get_db()
    session = next(session_generator)
    submit_github_issues(session, annotations)

    return "https://github.com/BTheDragonMaster/parasect/issues"


