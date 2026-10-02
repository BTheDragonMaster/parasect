# -*- coding: utf-8 -*-

"""Unit tests for parasect.core.retrain_models.model_needs_retraining, which
prepare_model() (parasect.core.helpers, on the `paras` model-preparation path)
consults before deciding whether to download or retrain a model.

Out of scope: retrain_model(), update_fingerprints(), update_metadata_file()'s
write side (covered via write_model_metadata_file in test_writers.py) -- these
require a populated PARASECT sqlite database and full training data and are
training-pipeline concerns, not part of running a `paras` prediction.
"""

import os
import unittest
from tempfile import TemporaryDirectory

import sklearn

from parasect.core.models import ModelType
from parasect.core.retrain_models import model_needs_retraining
from parasect.core.writers import write_model_metadata_file


class TestModelNeedsRetraining(unittest.TestCase):
    def _write_metadata(self, tmp_dir, model_to_version):
        path = os.path.join(tmp_dir, "metadata.txt")
        write_model_metadata_file(model_to_version, path)
        return path

    def test_returns_false_when_recorded_version_matches_installed_sklearn(self):
        # Compare against the real, currently-installed sklearn.__version__ rather
        # than a hardcoded version string, so this test is correct regardless of
        # which sklearn version is installed when it runs.
        with TemporaryDirectory() as tmp_dir:
            path = self._write_metadata(tmp_dir, {ModelType.PARAS: sklearn.__version__})
            self.assertFalse(model_needs_retraining(path, ModelType.PARAS))

    def test_returns_true_when_recorded_version_differs_from_installed_sklearn(self):
        with TemporaryDirectory() as tmp_dir:
            path = self._write_metadata(tmp_dir, {ModelType.PARAS: "0.0.1-not-a-real-version"})
            self.assertTrue(model_needs_retraining(path, ModelType.PARAS))

    def test_raises_for_composite_model_type(self):
        with TemporaryDirectory() as tmp_dir:
            path = self._write_metadata(tmp_dir, {ModelType.PARAS: sklearn.__version__})
            with self.assertRaises(ValueError):
                model_needs_retraining(path, ModelType.ANTISMASH_MODELS)

    def test_raises_keyerror_when_model_type_missing_from_metadata(self):
        with TemporaryDirectory() as tmp_dir:
            path = self._write_metadata(tmp_dir, {ModelType.PARAS: sklearn.__version__})
            with self.assertRaises(KeyError):
                model_needs_retraining(path, ModelType.PARASECT)


if __name__ == "__main__":
    unittest.main()
