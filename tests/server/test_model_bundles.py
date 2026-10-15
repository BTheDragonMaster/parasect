"""Exercise the webapp model-init and lazy-loader paths without Flask/Redis."""

import gzip
import hashlib
import io
import json
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import joblib
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "app" / "src" / "server"))
import fetch_models
from routes.model_loader import ModelSpec, MultiModelLoader
from parasect.core.model_bundle import validate_bundle, validate_fingerprints


FP = b"substrate_name\tsmiles\t11\t22\nexample\tC\t1\t0\n"


def release(model, fp_name):
    output = io.BytesIO()
    joblib.dump(SimpleNamespace(n_features_in_=512), output)
    content = {model.archive: gzip.compress(output.getvalue()), fp_name: FP}
    files = {name: {"key": name, "size": len(data), "checksum": "md5:" + hashlib.md5(data).hexdigest(),
                   "links": {"self": "https://example.test/" + name}}
             for name, data in content.items()}
    return files, content


@pytest.mark.parametrize("key,fp_name", [("parasect", "fingerprints.txt"),
                                       ("parasectBacterial", "fingerprints_bacterial.txt")])
def test_init_fetches_verifies_and_loads_matching_pair(tmp_path, key, fp_name):
    model = fetch_models.MODELS[key]
    files, content = release(model, fp_name)
    with patch("fetch_models.urlopen", side_effect=lambda request, **kw: io.BytesIO(
            content[request.full_url.rsplit("/", 1)[-1]])):
        entry = fetch_models.ensure(model, str(tmp_path), "23262276", files)
    manifest = {model.target: entry}
    fetch_models.verify([model], str(tmp_path), manifest)
    fetch_models.write_manifest(str(tmp_path), manifest)
    assert entry["fingerprints_file"] == fp_name
    assert entry["fingerprints_md5"] == hashlib.md5(FP).hexdigest()
    assert fetch_models.present_and_pinned(model, str(tmp_path), "23262276", manifest, False)
    assert not fetch_models.present_and_pinned(model, str(tmp_path), "another-release", manifest, False)

    loader = MultiModelLoader({key: ModelSpec("test", str(tmp_path / model.target),
                                            fingerprints_file=str(tmp_path / fp_name))})
    loaded = loader.get(key)
    assert loader.get(key) is loaded
    assert loaded.parasect_fingerprints_sha256_ == entry["fingerprints_sha256"]

    # Same byte count and feature width; only the column order changes.
    (tmp_path / fp_name).write_bytes(FP.replace(b"11\t22", b"22\t11"))
    assert not fetch_models.present_and_pinned(model, str(tmp_path), "23262276", manifest, False)
    with pytest.raises(ValueError, match="does not match"):
        validate_fingerprints(loaded, str(tmp_path / fp_name))
    # Disabling unpickling must not disable model/fingerprint pairing checks.
    with patch.dict("os.environ", {"PARAS_MODELS_VERIFY": "0"}), \
            pytest.raises(ValueError, match="fingerprint checksum mismatch"):
        fetch_models.verify([model], str(tmp_path), manifest)


def test_missing_sidecar_in_release_fails_before_replacing_existing_files(tmp_path):
    model = fetch_models.MODELS["parasect"]
    (tmp_path / model.target).write_bytes(b"old-model")
    (tmp_path / "fingerprints.txt").write_bytes(FP)
    files, _ = release(model, "fingerprints.txt")
    del files["fingerprints.txt"]
    with patch("fetch_models.download") as download, pytest.raises(SystemExit, match="missing fingerprints.txt"):
        fetch_models.ensure(model, str(tmp_path), "23262276", files)
    download.assert_not_called()
    assert (tmp_path / model.target).read_bytes() == b"old-model"
    assert (tmp_path / "fingerprints.txt").read_bytes() == FP


def test_corrupt_sidecar_download_preserves_previous_pair(tmp_path):
    model = fetch_models.MODELS["parasect"]
    files, content = release(model, "fingerprints.txt")
    content["fingerprints.txt"] = b"corrupt"
    (tmp_path / model.target).write_bytes(b"old-model")
    (tmp_path / "fingerprints.txt").write_bytes(FP)
    with patch("fetch_models.urlopen", side_effect=lambda request, **kw: io.BytesIO(
            content[request.full_url.rsplit("/", 1)[-1]])), \
            pytest.raises(SystemExit, match="Checksum mismatch"):
        fetch_models.ensure(model, str(tmp_path), "23262276", files)
    assert (tmp_path / model.target).read_bytes() == b"old-model"
    assert (tmp_path / "fingerprints.txt").read_bytes() == FP


def test_old_manifest_for_current_record_is_not_a_complete_bundle(tmp_path):
    model = fetch_models.MODELS["parasect"]
    (tmp_path / model.target).write_bytes(b"model")
    old_manifest = {model.target: {"record": "23262276", "size": 5, "archive_md5": "old"}}
    assert not fetch_models.present_and_pinned(model, str(tmp_path), "23262276", old_manifest, False)


def test_keep_existing_cannot_silently_keep_unpaired_parasect(tmp_path):
    model = fetch_models.MODELS["parasect"]
    (tmp_path / model.target).write_bytes(b"model")
    with pytest.raises(SystemExit, match="cannot keep an unverified PARASECT"):
        fetch_models.present_and_pinned(model, str(tmp_path), "23262276", {}, True)


def test_loader_rejects_unpaired_model_before_unpickling(tmp_path):
    model_path = tmp_path / "model.parasect"
    model_path.write_bytes(b"model")
    loader = MultiModelLoader({"test": ModelSpec("test", str(model_path),
                                               fingerprints_file=str(tmp_path / "fingerprints.txt"))})
    with patch.object(loader, "_default_load") as unpickle, pytest.raises(ValueError, match="Missing PARASECT"):
        loader.get("test")
    unpickle.assert_not_called()


def test_equal_size_model_replacement_is_rejected(tmp_path):
    model = fetch_models.MODELS["parasect"]
    files, content = release(model, "fingerprints.txt")
    with patch("fetch_models.urlopen", side_effect=lambda request, **kw: io.BytesIO(
            content[request.full_url.rsplit("/", 1)[-1]])):
        entry = fetch_models.ensure(model, str(tmp_path), "23262276", files)
    model_path = tmp_path / model.target
    data = model_path.read_bytes()
    model_path.write_bytes(bytes([data[0] ^ 1]) + data[1:])
    with pytest.raises(ValueError, match="model checksum mismatch"):
        validate_bundle(str(model_path), entry)
