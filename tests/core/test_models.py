# -*- coding: utf-8 -*-

"""Unit tests for parasect.core.models (ModelType, Model)."""

import os
import unittest
from tempfile import TemporaryDirectory

from joblib import load
from sklearn.ensemble import RandomForestClassifier

from parasect.core.models import Model, ModelType


class TestModelType(unittest.TestCase):
    """Tests for the ModelType IntFlag enum and its composite members."""

    def test_composite_flags_are_bitwise_unions(self):
        self.assertEqual(ModelType.ANTISMASH_MODELS, ModelType.PARAS_ALL_SUBSTRATES | ModelType.PARASECT_BACTERIAL)
        self.assertEqual(
            ModelType.ALL_MODELS,
            ModelType.PARAS | ModelType.PARASECT | ModelType.PARAS_ALL_SUBSTRATES | ModelType.PARASECT_BACTERIAL,
        )

    def test_single_model_membership_in_composites(self):
        self.assertIn(ModelType.PARAS, ModelType.ALL_MODELS)
        self.assertIn(ModelType.PARASECT, ModelType.ALL_MODELS)
        self.assertIn(ModelType.PARAS_ALL_SUBSTRATES, ModelType.ANTISMASH_MODELS)
        self.assertIn(ModelType.PARASECT_BACTERIAL, ModelType.ANTISMASH_MODELS)
        self.assertNotIn(ModelType.PARAS, ModelType.ANTISMASH_MODELS)

    def test_bit_count_distinguishes_single_from_composite_types(self):
        # parasect.core.retrain_models.model_needs_retraining() relies on
        # bit_count() == 1 to reject composite model types -- this pins that invariant.
        for single in (ModelType.PARAS, ModelType.PARASECT, ModelType.PARAS_ALL_SUBSTRATES, ModelType.PARASECT_BACTERIAL):
            self.assertEqual(single.bit_count(), 1)
        for composite in (ModelType.ANTISMASH_MODELS, ModelType.ALL_MODELS):
            self.assertGreater(composite.bit_count(), 1)


class TestModel(unittest.TestCase):
    """Tests for the Model dataclass, which wraps a fitted classifier for saving to disk."""

    @classmethod
    def setUpClass(cls):
        # A tiny but real fitted classifier (not a mock), so save()/load() round-trips
        # an actual RandomForestClassifier the way prepare_model() does in production.
        classifier = RandomForestClassifier(n_estimators=2, random_state=0)
        classifier.fit([[0, 0], [1, 1]], [0, 1])
        cls.classifier = classifier

    def test_save_creates_output_directory_and_dumps_loadable_model(self):
        model = Model(self.classifier, "test_model.gz")
        with TemporaryDirectory() as tmp_dir:
            out_dir = os.path.join(tmp_dir, "does_not_exist_yet")
            model.save(out_dir)

            out_path = os.path.join(out_dir, "test_model.gz")
            self.assertTrue(os.path.exists(out_path))

            loaded = load(out_path)
            self.assertEqual(loaded.predict([[0, 0]]).tolist(), self.classifier.predict([[0, 0]]).tolist())

    def test_save_does_not_error_when_output_directory_already_exists(self):
        model = Model(self.classifier, "test_model.gz")
        with TemporaryDirectory() as tmp_dir:
            # tmp_dir already exists; save() must not try to mkdir it again.
            model.save(tmp_dir)
            self.assertTrue(os.path.exists(os.path.join(tmp_dir, "test_model.gz")))


if __name__ == "__main__":
    unittest.main()
