# -*- coding: utf-8 -*-

"""API for PARSECT."""

from flask import Response, jsonify
from routes.app import app
from routes.retrieval import blueprint_retrieve
from routes.submit import blueprint_submit_raw, blueprint_submit_quick
from routes.data_annotation import blueprint_annotate_data
from routes.annotation_editor import blueprint_check_smiles, blueprint_check_substrate_name, blueprint_get_substrates, \
    blueprint_submit_annotations, blueprint_check_protein_name, blueprint_check_domain_name
from routes.database import engine

from parasect.version import get_version

app.register_blueprint(blueprint_retrieve)
app.register_blueprint(blueprint_submit_raw)
app.register_blueprint(blueprint_submit_quick)
app.register_blueprint(blueprint_annotate_data)
app.register_blueprint(blueprint_check_smiles)
app.register_blueprint(blueprint_check_domain_name)
app.register_blueprint(blueprint_check_substrate_name)
app.register_blueprint(blueprint_get_substrates)
app.register_blueprint(blueprint_submit_annotations)
app.register_blueprint(blueprint_check_protein_name)


@app.errorhandler(404)
def not_found(error) -> Response:
    return app.send_static_file("index.html")


@app.route("/")
def index() -> Response:
    return app.send_static_file("index.html")


# api endpoint for fetching version
@app.route("/api/version")
def version() -> Response:
    return jsonify({"version": get_version()})
