# -*- coding: utf-8 -*-

"""Module for defining the Flask app."""

from flask import Flask

app = Flask(__name__)
app.config["JOB_RESULTS"] = dict()
