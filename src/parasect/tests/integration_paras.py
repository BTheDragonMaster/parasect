import unittest

import os
from tempfile import TemporaryDirectory
import argparse
import logging
from joblib import load
from shutil import copy

from parasect.core.constants import SEPARATOR_1, SEPARATOR_2, SEPARATOR_3
from parasect.api import run_paras
from parasect.core.helpers import download_and_unpack_or_fetch
from parasect.core.writers import write_fasta_file, write_results
from parasect.core.retrain_models import retrain_model, model_needs_retraining, update_metadata_file
from parasect.core.models import ModelType
from parasect.core.constants import MODEL_METADATA_FILE



class TestParas(unittest.TestCase):
    def setUp(self):
        self.data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
        self.temp_dir = TemporaryDirectory()

    def testDptA(self):
        input_file = os.path.join(self.data_dir, "dptA.fasta")
        with open(input_file) as f:
            domain_data = f.read()
            run_paras(domain_data, "fasta", self.temp_dir, )

class TestParasect(unittest.TestCase):
    pass