"""Keep PARASECT models paired with their ordered substrate fingerprints."""

import hashlib
import os


FINGERPRINT_FILES = {
    "model.parasect": "fingerprints.txt",
    "bacterial_model.parasect": "fingerprints_bacterial.txt",
}


def fingerprint_path(model_path: str) -> str:
    """Return the fingerprint sidecar path for a PARASECT model."""
    name = os.path.basename(model_path).removesuffix(".gz")
    return os.path.join(os.path.dirname(model_path), FINGERPRINT_FILES[name])


def file_sha256(path: str) -> str:
    """Hash a file without loading the model into memory."""
    digest = hashlib.sha256()
    with open(path, "rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bundle_checksums(model_path: str) -> dict:
    """Record both files after downloading or training them together."""
    fingerprints = fingerprint_path(model_path)
    return {
        "model_sha256": file_sha256(model_path),
        "fingerprints_file": os.path.basename(fingerprints),
        "fingerprints_sha256": file_sha256(fingerprints),
    }


def validate_bundle(model_path: str, metadata: dict) -> str:
    """Reject missing, replaced, or mismatched model/fingerprint files."""
    fingerprints = fingerprint_path(model_path)
    if not metadata or metadata.get("fingerprints_file") != os.path.basename(fingerprints):
        raise ValueError("Missing PARASECT model/fingerprint metadata; prepare the model again.")
    if file_sha256(model_path) != metadata.get("model_sha256"):
        raise ValueError("PARASECT model checksum mismatch; prepare the model again.")
    if file_sha256(fingerprints) != metadata.get("fingerprints_sha256"):
        raise ValueError("PARASECT fingerprint checksum mismatch; prepare the model again.")
    return fingerprints


def validate_fingerprints(model, fingerprints_file: str, expected_sha256: str = None) -> None:
    """Check schema width and, when available, the exact training fingerprint file.

    Width alone cannot catch reordered columns (issue #100). Managed loaders
    bind the verified file's checksum to the in-memory model so predictions also
    reject a sidecar that changes after the model has been cached.
    """
    expected = expected_sha256 or getattr(model, "parasect_fingerprints_sha256_", None)
    if expected and file_sha256(fingerprints_file) != expected:
        raise ValueError("Fingerprint file does not match the loaded PARASECT model.")
    with open(fingerprints_file) as source:
        header = source.readline().rstrip("\r\n").split("\t")
    if header[:2] != ["substrate_name", "smiles"] or len(header) <= 2:
        raise ValueError("Invalid PARASECT fingerprint header.")
    hashes = [int(value) for value in header[2:]]
    if len(set(hashes)) != len(hashes):
        raise ValueError("Duplicate columns in PARASECT fingerprint file.")
    # Extended domain signatures contain 34 residues, with 15 properties each.
    if model.n_features_in_ != 34 * 15 + len(hashes):
        raise ValueError("Fingerprint feature count does not match the PARASECT model.")
    if expected:
        model.parasect_fingerprints_sha256_ = expected
