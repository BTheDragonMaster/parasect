# -*- coding: utf-8 -*-

"""API for PARSECT."""

from flask import Response
from routes.app import app
from routes.retrieval import blueprint_retrieve
from routes.submit_raw import blueprint_submit_raw
from routes.submit_signature import blueprint_submit_signature

from parasect.version import get_version

app.register_blueprint(blueprint_retrieve)
app.register_blueprint(blueprint_submit_raw)
app.register_blueprint(blueprint_submit_signature)


@app.errorhandler(404)
def not_found(error) -> Response:
    return app.send_static_file("index.html")


@app.route("/")
def index() -> Response:
    return app.send_static_file("index.html")


# api endpoint for fetching version
@app.route("/api/version")
def version() -> Response:
    return {"version": get_version()}


def main() -> None:
    app.run(host="localhost", port=4000, debug=True)


if __name__ == "__main__":
    main()
