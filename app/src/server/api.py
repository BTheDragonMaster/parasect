# -*- coding: utf-8 -*-

"""API for PARSECT."""

from flask import Flask, Response
from routes.submit import blueprint_submit_raw

from parasect.version import get_version

app = Flask(__name__)
app.register_blueprint(blueprint_submit_raw)


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
