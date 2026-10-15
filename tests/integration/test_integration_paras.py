"""
Integration tests for PARAS/PARASECT models. Only use with pytest.
"""

import unittest
import pytest

import os
from tempfile import TemporaryDirectory
from joblib import load

from parasect.api import run_paras, run_parasect, Result
from parasect.core.helpers import prepare_model, load_parasect_bundle
from parasect.core.writers import write_results
from parasect.core.models import ModelType


class BaseModelTest(unittest.TestCase):
    """Base testing class for PARAS/PARASECT models. Not tested itself."""
    __test__ = False
    MODEL_TYPE = None

    @classmethod
    def setUpClass(cls):
        """Prepares directories and loads the model"""
        if cls.MODEL_TYPE is None:
            raise unittest.SkipTest(f"{cls.__name__} is a base class, not run directly")

        cls.data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

        cls.temp_dir = TemporaryDirectory()
        cls.addClassCleanup(cls.temp_dir.cleanup)
        cls.model_dir = TemporaryDirectory()
        cls.addClassCleanup(cls.model_dir.cleanup)
        cls.output_dir = TemporaryDirectory()
        cls.addClassCleanup(cls.output_dir.cleanup)

        model_path = prepare_model(cls.MODEL_TYPE, cls.model_dir.name)
        if cls.MODEL_TYPE in (ModelType.PARASECT, ModelType.PARASECT_BACTERIAL):
            cls.model, cls.fingerprints_file = load_parasect_bundle(model_path)
        else:
            cls.model = load(model_path)

    def _read_input(self, file_name: str) -> str:
        input_file = os.path.join(self.data_dir, file_name)
        with open(input_file) as f:
            return f.read()

    def _write_results(self, results: list[Result], n_predictions: int, **write_kwargs) -> None:
        write_results(results,
                      self.output_dir.name,
                      n_predictions,
                      self.MODEL_TYPE,
                      **write_kwargs)

    def _run_and_save(self, job_name: str, file_type: str = "fasta") -> list[Result]:
        domain_data = self._read_input(f"{job_name}.{file_type}")
        results = self._run_model(domain_data)
        self._write_results(results, 3,
                            job_name=job_name,
                            save_signatures=True,
                            save_extended_signatures=True,
                            save_domain_sequences=True)
        return results

    def _assert_output_exists(self, job_name: str,
                              suffixes=("signatures", "extended_signatures", "sequences")):
        for suffix in suffixes:
            path = os.path.join(self.output_dir.name, f"{job_name}_{suffix}.fasta")
            self.assertTrue(os.path.exists(path), f"Missing expected output: {path}")

    def _run_model(self, domain_data: str):
        raise NotImplementedError()


class ParasBaseTest(BaseModelTest):
    __test__ = False
    def _run_model(self, domain_data: str) -> list[Result]:
        results = run_paras(domain_data,
                            "fasta",
                            self.temp_dir.name,
                            self.model)

        return results


class ParasectBaseTest(BaseModelTest):
    __test__ = False
    def _run_model(self, domain_data: str):
        results = run_parasect(domain_data,
                               "fasta",
                               self.temp_dir.name,
                               self.model,
                               bacterial_only=self.MODEL_TYPE == ModelType.PARASECT_BACTERIAL,
                               fingerprints_file=self.fingerprints_file)

        return results


@pytest.mark.integration
class TestParas(ParasBaseTest):
    __test__ = True
    MODEL_TYPE = ModelType.PARAS

    def testDptA(self):
        job_name = "dptA"
        results = self._run_and_save(job_name)
        self.assertEqual(len(results), 5)
        self._assert_output_exists(job_name)

    def testFragmentedHit(self):
        job_name = "fragmented_hit"
        results = self._run_and_save(job_name)
        self.assertEqual(len(results), 1)
        self._assert_output_exists(job_name)

    def testPoorQualityDomain(self):
        job_name = "poor_quality_domain"
        results = self._run_and_save(job_name)
        self.assertEqual(len(results), 1)
        self._assert_output_exists(job_name)

    def testAoxDomain(self):
        job_name = "a_ox"
        results = self._run_and_save(job_name)
        self.assertEqual(len(results), 1)
        self._assert_output_exists(job_name)

    def testOxDomain(self):
        job_name = "ox"
        with self.assertRaises(Exception) as error:
            self._run_and_save(job_name)
        self.assertEqual(str(error.exception), "no adenylation domains found")


@pytest.mark.integration
class TestParasAllSubstrates(ParasBaseTest):
    __test__ = True
    MODEL_TYPE = ModelType.PARAS_ALL_SUBSTRATES

    def testDptA(self):
        job_name = "dptA"
        results = self._run_and_save(job_name)
        self.assertEqual(len(results), 5)
        self._assert_output_exists(job_name)

    def testFragmentedHit(self):
        job_name = "fragmented_hit"
        results = self._run_and_save(job_name)
        self.assertEqual(len(results), 1)
        self._assert_output_exists(job_name)

    def testPoorQualityDomain(self):
        job_name = "poor_quality_domain"
        results = self._run_and_save(job_name)
        self.assertEqual(len(results), 1)
        self._assert_output_exists(job_name)

    def testAoxDomain(self):
        job_name = "a_ox"
        results = self._run_and_save(job_name)
        self.assertEqual(len(results), 1)
        self._assert_output_exists(job_name)

    def testOxDomain(self):
        job_name = "ox"
        with self.assertRaises(Exception) as error:
            self._run_and_save(job_name)
        self.assertEqual(str(error.exception), "no adenylation domains found")

@pytest.mark.integration
class TestParasect(ParasectBaseTest):
    __test__ = True
    MODEL_TYPE = ModelType.PARASECT

    def test_dpta_standard_and_custom_substrates_agree(self):
        results = self._run_model(self._read_input("dptA.fasta"))
        self.assertEqual(len(results), 5)
        self._assert_custom_predictions_match(results)

    def _assert_custom_predictions_match(self, results):
        reference = results[0]
        custom = run_parasect(
            self._read_input("dptA.fasta"), "fasta", self.temp_dir.name, self.model,
            custom_substrate_names=reference.prediction_labels,
            custom_substrate_smiles=[row["substrate_smiles"] for row in reference.to_json()["predictions"]],
            only_custom=True,
            bacterial_only=self.MODEL_TYPE == ModelType.PARASECT_BACTERIAL,
            fingerprints_file=self.fingerprints_file,
        )
        for standard_result, custom_result in zip(results, custom):
            self.assertEqual(standard_result.predictions, custom_result.predictions)

@pytest.mark.integration
class TestParasectBacterial(TestParasect):
    MODEL_TYPE = ModelType.PARASECT_BACTERIAL
