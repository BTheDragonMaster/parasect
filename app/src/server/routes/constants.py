# -*- coding: utf-8 -*-

"""Constants used throughout the server package."""

import os

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
INCLUDED_SUBSTRATES_FILE = os.path.join(DATA_DIR, "included_substrates.txt")
FINGERPRINTS_FILE = os.path.join(DATA_DIR, "fingerprints.txt")

MODEL_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "models")
TEMP_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "temp")
