# -*- coding: utf-8 -*-

"""Constants used throughout the server package."""

import os

MODEL_DIR_LOCAL = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))), "models")
MODEL_DIR = os.getenv("MODEL_DIR", MODEL_DIR_LOCAL)
TEMP_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "temp")
