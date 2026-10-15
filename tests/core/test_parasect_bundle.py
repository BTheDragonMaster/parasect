"""Regressions for issue #100: feature order belongs to the model release."""

import hashlib
import io
import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest

from parasect.api import run_parasect
from parasect.core.helpers import prepare_model, load_parasect_bundle
from parasect.core.model_bundle import (
    bundle_checksums, fingerprint_path, validate_bundle, validate_fingerprints,
)
from parasect.core.models import Model, ModelType
from parasect.core.retrain_models import retrain_model


FP = b"substrate_name\tsmiles\t11\t22\nexample\tC\t1\t0\n"


def release_files(model_name, fp_name, model_bytes=b"model", fp_bytes=FP):
    return {name: {"key": name, "checksum": "md5:" + hashlib.md5(content).hexdigest(),
                   "size": len(content), "links": {"self": "https://example.test/" + name}}
            for name, content in [(model_name, model_bytes), (fp_name, fp_bytes)]}


@pytest.mark.parametrize("model_type,name,fp_name", [
    (ModelType.PARASECT, "model.parasect.gz", "fingerprints.txt"),
    (ModelType.PARASECT_BACTERIAL, "bacterial_model.parasect.gz", "fingerprints_bacterial.txt"),
])
def test_cli_replaces_legacy_cache_and_reuses_verified_pair(tmp_path, model_type, name, fp_name):
    model_path = tmp_path / name
    model_path.write_bytes(b"old-model")
    (tmp_path / fp_name).write_bytes(b"old-fingerprints")
    files = release_files(name, fp_name)

    def response(request, **kwargs):
        url = request.full_url
        if "/api/records/" in url:
            return io.BytesIO(json.dumps({"files": list(files.values())}).encode())
        return io.BytesIO(FP if url.endswith(fp_name) else b"model")

    with patch("parasect.core.helpers.model_needs_retraining", return_value=False), \
            patch("parasect.core.helpers.urlopen", side_effect=response) as network:
        assert prepare_model(model_type, str(tmp_path)) == str(model_path)
        assert model_path.read_bytes() == b"model"
        assert (tmp_path / fp_name).read_bytes() == FP
        assert network.call_count == 3
        network.reset_mock()
        prepare_model(model_type, str(tmp_path))
        network.assert_not_called()
        # A release change must invalidate the pair even when filenames match.
        with patch("parasect.core.helpers.ZENODO_RECORD", "next-release"):
            prepare_model(model_type, str(tmp_path))
        assert network.call_count == 3


def test_failed_fingerprint_download_preserves_previous_cli_pair(tmp_path):
    model_path = tmp_path / "model.parasect.gz"
    model_path.write_bytes(b"old-model")
    fp_path = Path(fingerprint_path(str(model_path)))
    fp_path.write_bytes(FP)
    files = release_files(model_path.name, fp_path.name)

    def response(request, **kwargs):
        if "/api/records/" in request.full_url:
            return io.BytesIO(json.dumps({"files": list(files.values())}).encode())
        return io.BytesIO(b"corrupt" if request.full_url.endswith(fp_path.name) else b"model")

    with patch("parasect.core.helpers.model_needs_retraining", return_value=False), \
            patch("parasect.core.helpers.urlopen", side_effect=response), \
            pytest.raises(ValueError, match="Checksum mismatch"):
        prepare_model(ModelType.PARASECT, str(tmp_path))
    assert model_path.read_bytes() == b"old-model"
    assert fp_path.read_bytes() == FP
    assert not Path(str(model_path) + ".json").exists()


@pytest.mark.parametrize("model_type,filename", [
    (ModelType.PARASECT, "model.parasect.gz"),
    (ModelType.PARASECT_BACTERIAL, "bacterial_model.parasect.gz"),
])
def test_cli_retraining_keeps_fingerprints_with_the_model(tmp_path, model_type, filename):
    def train(selected_type, fingerprints_out_file):
        assert selected_type == model_type
        assert Path(fingerprints_out_file).is_relative_to(tmp_path)
        Path(fingerprints_out_file).write_bytes(FP)
        return Model(SimpleNamespace(n_features_in_=512), filename)

    with patch("parasect.core.helpers.model_needs_retraining", return_value=True), \
            patch("parasect.core.helpers.retrain_model", side_effect=train) as training, \
            patch("parasect.core.helpers.urlopen") as network:
        model_path = prepare_model(model_type, str(tmp_path))
        model, fp_path = load_parasect_bundle(model_path)
        assert model.n_features_in_ == 512
        assert Path(fp_path).read_bytes() == FP
        prepare_model(model_type, str(tmp_path))
        training.assert_called_once()
        network.assert_not_called()


@pytest.mark.parametrize("model_type", [ModelType.PARASECT, ModelType.PARASECT_BACTERIAL])
def test_retraining_requires_explicit_output_instead_of_overwriting_package(model_type):
    with pytest.raises(ValueError, match="fingerprint output path"):
        retrain_model(model_type)


def test_equal_width_reordered_fingerprints_are_rejected(tmp_path):
    model_path = tmp_path / "model.parasect"
    model_path.write_bytes(b"model")
    fp_path = Path(fingerprint_path(str(model_path)))
    fp_path.write_bytes(FP)
    metadata = bundle_checksums(str(model_path))
    model = SimpleNamespace(n_features_in_=512)
    validate_fingerprints(model, str(fp_path), metadata["fingerprints_sha256"])
    fp_path.write_bytes(FP.replace(b"11\t22", b"22\t11"))
    with pytest.raises(ValueError, match="fingerprint checksum mismatch"):
        validate_bundle(str(model_path), metadata)
    with pytest.raises(ValueError, match="does not match"):
        validate_fingerprints(model, str(fp_path))


@pytest.mark.parametrize("bacterial_only", [False, True])
@pytest.mark.parametrize("only_custom", [False, True])
def test_prediction_uses_explicit_file_for_standard_and_custom(tmp_path, bacterial_only, only_custom):
    # Reversed feature order is deliberate. Both routes must obey this file,
    # rather than either packaged fingerprint constant.
    fp_path = tmp_path / "custom-order.tsv"
    fp_path.write_text("substrate_name\tsmiles\t22\t11\nexample\tC\t0\t1\n")
    seen = []

    class RecordingModel:
        n_features_in_ = 512

        def predict_proba(self, vectors):
            seen.extend(vectors)
            return np.array([[0.2, 0.8]] * len(vectors))

    with patch("parasect.api.parse_substrate_list", return_value=["example"]), \
            patch("parasect.api.get_domains", return_value=[SimpleNamespace(extended_signature="A" * 34)]), \
            patch("parasect.core.parsing.smiles_to_fingerprint", return_value={11}):
        results = run_parasect(
            ">test\nAAAA", "fasta", str(tmp_path), RecordingModel(),
            custom_substrate_names=["custom"], custom_substrate_smiles=["C"],
            only_custom=only_custom, bacterial_only=bacterial_only, fingerprints_file=str(fp_path),
        )
    assert len(seen) == (1 if only_custom else 2)
    assert all(vector[-2:] == [0, 1] for vector in seen)
    assert results[0].prediction_labels == (["custom"] if only_custom else ["example", "custom"])


def test_api_requires_explicit_fingerprints(tmp_path):
    with pytest.raises(TypeError, match="fingerprints_file"):
        run_parasect(">test\nAAAA", "fasta", str(tmp_path), SimpleNamespace(n_features_in_=512))


def test_wrong_feature_count_fails_before_domain_extraction(tmp_path):
    fp_path = tmp_path / "fingerprints.txt"
    fp_path.write_bytes(FP)
    with patch("parasect.api.get_domains") as domains, pytest.raises(ValueError, match="feature count"):
        run_parasect(">test\nAAAA", "fasta", str(tmp_path), SimpleNamespace(n_features_in_=513),
                     fingerprints_file=str(fp_path))
    domains.assert_not_called()
