# -*- coding: utf-8 -*-

"""Module for defining the Flask app."""

import logging
import os
import time

from flask import Flask, request


logging.basicConfig(level=os.getenv("LOG_LEVEL", "INFO"))


app = Flask(__name__)


@app.before_request
def _t0():
    request._t0 = time.time()


@app.after_request
def _log(resp):
    dt = (time.time() - getattr(request, "_t0", time.time())) * 1000
    app.logger.info("method=%s path=%s status=%s dt_ms=%.1f", request.method, request.path, resp.status_code, dt)
    return resp


@app.get("/health")
def health():
    from .job_store import ping as redis_ping

    if not redis_ping():
        return {"ok": False, "redis": False}, 503
    return {"ok": True, "redis": True}, 200


app.config["ENV"] = os.getenv("FLASK_ENV", "production")  # defaults to "production"
app.config["DEBUG"] = app.config["ENV"] == "development"
print("starting app in environment:", app.config["ENV"])
print("Debug mode is:", app.debug)


if app.config["ENV"] == "production":
    print("production environment detected")
elif app.config["ENV"] == "development":
    print("development environment detected")
else:
    print(f"unknown environment: {app.config['ENV']}")
