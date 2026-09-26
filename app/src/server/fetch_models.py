# -*- coding: utf-8 -*-

"""Fetch the PARAS/PARASECT model files this server needs into MODEL_DIR.

Run as a one-shot init step before the server starts (see the ``paras-models``
service in docker-compose.yml), or by hand for local development:

    python fetch_models.py

Models are published on Zenodo as gzipped pickles (~214 MB for the full set),
but the server needs them uncompressed: routes/submit.py loads them with
mmap=True so gunicorn workers share pages through the OS page cache, and
joblib silently ignores mmap_mode on a compressed file. Expanding them is
what turns that download into ~4 GB on disk.

Why not parasect.core.helpers.prepare_model: that consults
model_metadata.txt, which records the exact scikit-learn version each model
was pickled with. When the running interpreter has any other version it takes
the retrain branch and trains four random forests from scratch. Reasonable
for a CLI run on a workstation; completely wrong for a container coming up
behind a health check!

Environment:
    MODEL_DIR       where to put the files (default: app/models beside this
                    checkout; docker-compose sets /app/models).
    PARAS_MODELS    comma-separated model keys to fetch, using the same names
                    the API and UI use (see MODELS below), or ``all``.
                    Default: all of them.
    ZENODO_RECORD   Zenodo record id to pin to (default below).
    PARAS_MODELS_KEEP_EXISTING
                    if set to 1/true/yes, leave already-present files alone
                    even when they didn't come from the pinned record. For
                    when you are deliberately running your own weights.
    PARAS_MODELS_VERIFY
                    set to 0/false/no to skip the post-fetch load check.
"""

from __future__ import annotations

import errno
import gzip
import hashlib
import json
import logging
import os
import shutil
import sys
from typing import NamedTuple
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


# Zenodo record 18682178 ("v5", Feb 2026): the models re-pickled under
# scikit-learn 1.8.0.
#
# The same record parasect/core/helpers.py pins for the CLI (ZENODO_RECORD
# there). Kept as a separate constant so a deployment can be re-pinned without
# touching the core package, but the two should normally agree.
#
# It replaced record 17224548, whose models were pickled under sklearn 1.2.0:
# sklearn grew a field in the decision-tree node dtype in 1.3, so those fail to
# load on anything newer with "node array from the pickle has an incompatible
# dtype". This is an error, not a warning. Move this together with the scikit-learn
# pin in server-requirements.txt; they are one unit.
#
# Pinned rather than tracking Zenodo's "latest version" pointer: a floating
# pointer would let someone else's upload silently change what this deployment
# predicts, with no commit and no way to reproduce a result a user got last
# week. Bumping this should be a visible one-line change.
DEFAULT_ZENODO_RECORD = "18682178"

ZENODO_API = "https://zenodo.org/api/records/{record}"

USER_AGENT = "parasect-webapp model fetcher"

# Stream in chunks: the largest archive is ~132 MB compressed but expands past
# 3 GB, and neither end of that belongs in memory.
CHUNK = 1024 * 1024

# Written next to the models so a later run can tell which record a file on
# disk came from. Without it, "pinning" would be decorative: every file would
# look present and a changed pin would never take effect.
MANIFEST_NAME = ".fetch_manifest.json"


class ModelFile(NamedTuple):
    """One fetchable model: its key in the API, and its names on each side."""

    key: str          # what the client sends as `selectedModel`
    archive: str      # file name in the Zenodo record
    target: str       # file name routes/submit.py opens
    label: str        # for logs
    expanded: int     # approx. uncompressed bytes, for the disk-space check


# Keyed by the strings the submit endpoint accepts, so pinning the set of models
# here and the set offered in the UI are obviously the same list. See the
# model registry at the top of routes/submit.py.
#
# `expanded` is measured against the default record; it only feeds the
# pre-flight free-space warning, so it being a little stale is harmless.
MODELS: dict[str, ModelFile] = {
    m.key: m
    for m in (
        ModelFile("parasAllSubstrates", "all_substrates_model.paras.gz",
                  "all_substrates_model.paras", "PARAS (all substrates)", 3_473_913_870),
        ModelFile("parasCommonSubstrates", "model.paras.gz",
                  "model.paras", "PARAS (common substrates)", 321_878_398),
        ModelFile("parasect", "model.parasect.gz",
                  "model.parasect", "PARASECT", 98_986_566),
        ModelFile("parasectBacterial", "bacterial_model.parasect.gz",
                  "bacterial_model.parasect", "PARASECT (bacterial)", 74_848_646),
    )
}


def env_flag(name: str) -> bool:
    return (os.getenv(name) or "").strip().lower() in {"1", "true", "yes", "on"}


def env_flag_default_true(name: str) -> bool:
    """Like function env_flag, but unset means on."""
    raw = (os.getenv(name) or "").strip().lower()
    return raw not in {"0", "false", "no", "off"}


def default_model_dir() -> str:
    """``app/models`` relative to this checkout.

    Only for local runs; the container always sets MODEL_DIR explicitly,
    because in the image this file sits at /app/fetch_models.py and the walk
    below would land somewhere meaningless.
    """
    server_dir = os.path.dirname(os.path.abspath(__file__))
    app_dir = os.path.dirname(os.path.dirname(server_dir))
    return os.path.join(app_dir, "models")


def selected_models(raw: str) -> list[ModelFile]:
    """Resolve the PARAS_MODELS setting into model entries.

    An unknown key is an error rather than a warning: a typo would otherwise
    bring the server up looking healthy, and only fail at the first prediction
    request that happened to want the model nobody downloaded.
    """
    raw = (raw or "").strip()
    if not raw or raw.lower() == "all":
        return list(MODELS.values())

    wanted: list[ModelFile] = []
    unknown: list[str] = []
    for name in (part.strip() for part in raw.split(",")):
        if not name:
            continue
        if name in MODELS:
            if MODELS[name] not in wanted:
                wanted.append(MODELS[name])
        else:
            unknown.append(name)

    if unknown:
        raise SystemExit(
            f"PARAS_MODELS lists unknown model(s): {', '.join(unknown)}. "
            f"Valid keys: {', '.join(MODELS)} (or 'all')."
        )
    if not wanted:
        raise SystemExit("PARAS_MODELS is set but selects no models.")
    return wanted


def read_manifest(model_dir: str) -> dict[str, dict]:
    try:
        with open(os.path.join(model_dir, MANIFEST_NAME)) as fh:
            data = json.load(fh)
        return data.get("files", {}) if isinstance(data, dict) else {}
    except (OSError, ValueError):
        # Absent or unreadable is the normal first-run case, and a corrupt one
        # should just mean "re-verify everything", never a crash.
        return {}


def write_manifest(model_dir: str, files: dict[str, dict]) -> None:
    path = os.path.join(model_dir, MANIFEST_NAME)
    partial = f"{path}.part"
    with open(partial, "w") as fh:
        json.dump({"files": files}, fh, indent=2, sort_keys=True)
    os.replace(partial, path)


def present_and_pinned(model: ModelFile, model_dir: str, record: str,
                       manifest: dict[str, dict], keep_existing: bool) -> bool:
    """Decide whether model on disk can be left as it is.

    Three cases, because "the file exists" is not the same as "the file is the
    one we pinned":
      - we fetched it from this record  -> keep it
      - we fetched it from another one  -> the pin moved, replace it
      - we never fetched it             -> unverifiable, so replace it unless
                                           the operator says it's theirs
    """
    target = os.path.join(model_dir, model.target)
    if not (os.path.exists(target) and os.path.getsize(target) > 0):
        return False

    entry = manifest.get(model.target)
    if entry and entry.get("record") == record and entry.get("size") == os.path.getsize(target):
        logging.info("%s: already present from record %s", model.label, record)
        return True

    if entry:
        logging.info("%s: on disk from record %s, pinned to %s, replacing",
                     model.label, entry.get("record", "?"), record)
        return False

    if keep_existing:
        logging.warning(
            "%s: %s is already there but wasn't fetched by this script, so it "
            "can't be checked against record %s. Keeping it "
            "(PARAS_MODELS_KEEP_EXISTING is set).", model.label, model.target, record)
        return True

    logging.warning(
        "%s: %s is already there but wasn't fetched by this script, so it "
        "can't be checked against record %s, so replacing it. Set "
        "PARAS_MODELS_KEEP_EXISTING=1 to keep your own copy instead.",
        model.label, model.target, record)
    return False


def fetch_record(record: str) -> dict[str, dict]:
    """Return the record's files, keyed by file name.

    Checksums come from the record itself rather than being hard-coded, so
    re-pinning ZENODO_RECORD can't leave us verifying against the old hashes.
    """
    url = ZENODO_API.format(record=record)
    logging.info("Reading Zenodo record %s", record)
    try:
        with urlopen(Request(url, headers={"User-Agent": USER_AGENT}), timeout=60) as resp:
            payload = json.load(resp)
    except (HTTPError, URLError) as e:
        raise SystemExit(f"Could not read Zenodo record {record}: {e}")

    files = {f["key"]: f for f in payload.get("files", [])}
    if not files:
        raise SystemExit(f"Zenodo record {record} lists no files.")
    return files


def download(url: str, dest: str, expected_md5: str, size: int) -> None:
    """Stream ``url`` to ``dest``, failing if the md5 doesn't match.

    Hashes while writing rather than re-reading the finished file, and leaves
    nothing behind on a mismatch: a truncated archive would otherwise expand
    into a model that only fails much later, inside a prediction.
    """
    digest = hashlib.md5()
    done = 0
    next_report = 10

    try:
        with urlopen(Request(url, headers={"User-Agent": USER_AGENT}), timeout=60) as resp, \
                open(dest, "wb") as out:
            while True:
                chunk = resp.read(CHUNK)
                if not chunk:
                    break
                out.write(chunk)
                digest.update(chunk)
                done += len(chunk)
                if size:
                    pct = done * 100 / size
                    if pct >= next_report:
                        logging.info("    %d%% (%.1f/%.1f MB)", pct, done / 1e6, size / 1e6)
                        next_report += 10
    except (HTTPError, URLError, OSError) as e:
        if os.path.exists(dest):
            os.remove(dest)
        raise SystemExit(f"Download of {os.path.basename(dest)} failed: {e}")

    actual = digest.hexdigest()
    if actual != expected_md5:
        os.remove(dest)
        raise SystemExit(
            f"Checksum mismatch for {os.path.basename(dest)}: "
            f"expected md5 {expected_md5}, got {actual}. Refusing to use it."
        )


def expand(archive: str, target: str) -> None:
    """Gunzip ``archive`` into ``target``, atomically.

    The expansion is written beside the target and renamed into place only once
    complete, so it also replaces any previous file in one step. Without that, a
    container killed partway through (past 3 GB for the largest model, so the
    window is not small) would leave a truncated file that looks present to the
    next run and then blows up inside a prediction.
    """
    partial = f"{target}.part"
    try:
        with gzip.open(archive, "rb") as src, open(partial, "wb") as dst:
            shutil.copyfileobj(src, dst, CHUNK)
        os.replace(partial, target)
    except BaseException:
        if os.path.exists(partial):
            os.remove(partial)
        raise


def ensure(model: ModelFile, model_dir: str, record: str,
           record_files: dict[str, dict]) -> dict:
    """Download and expand one model. Returns its manifest entry."""
    entry = record_files.get(model.archive)
    if entry is None:
        raise SystemExit(
            f"{model.label}: the pinned Zenodo record has no file named "
            f"'{model.archive}'. Available: {', '.join(sorted(record_files))}."
        )

    expected_md5 = entry["checksum"].split(":", 1)[-1]  # Zenodo reports "md5:<hex>"
    archive = os.path.join(model_dir, f"{model.archive}.part")
    target = os.path.join(model_dir, model.target)

    logging.info("%s: downloading %s (%.1f MB)", model.label, model.archive, entry["size"] / 1e6)
    download(entry["links"]["self"], archive, expected_md5, entry.get("size", 0))

    logging.info("%s: expanding to %s", model.label, model.target)
    try:
        expand(archive, target)
    finally:
        # The compressed copy has served its purpose; keeping it is how a model
        # directory ends up holding both a .gz and its expansion, which is only
        # ever confusing.
        if os.path.exists(archive):
            os.remove(archive)

    size = os.path.getsize(target)
    logging.info("%s: ready (%.2f GB)", model.label, size / 1e9)
    return {"record": record, "archive_md5": expected_md5, "size": size}


def verify(models: list[ModelFile], model_dir: str, manifest: dict[str, dict]) -> None:
    """Actually load every requested model, and fail loudly if one won't.

    This is the check whose absence let a real breakage sit unnoticed: a random
    forest pickled under one scikit-learn and read under another can fail
    outright ("node array from the pickle has an incompatible dtype"), and
    because MultiModelLoader loads lazily, nothing noticed until a user
    submitted a prediction, by which time the server was long since up and
    reporting healthy.

    Run on every init rather than only after a download, because the pairing
    that breaks is (model file, scikit-learn version) and either side can move:
    a rebuild that resolves a different scikit-learn leaves the files untouched
    and still breaks them.

    Loading with mmap keeps this cheap in memory; it costs a few seconds, in a
    one-shot container whose whole job is to fail before the server starts.
    """
    if not env_flag_default_true("PARAS_MODELS_VERIFY"):
        logging.warning("PARAS_MODELS_VERIFY is off, and not checking that the "
                        "models actually load.")
        return

    try:
        import joblib
        import sklearn
    except ImportError as e:
        logging.warning("Skipping load check (%s). The server will be the first "
                        "thing to find out if these models don't load.", e)
        return

    logging.info("Checking the models load under scikit-learn %s", sklearn.__version__)
    broken: list[str] = []
    for model in models:
        target = os.path.join(model_dir, model.target)
        try:
            try:
                joblib.load(target, mmap_mode="r")
            except OSError as e:
                # Mirror MultiModelLoader._default_load in routes/model_loader.py:
                # mmap'ing a random forest costs a descriptor per array, so a low
                # RLIMIT_NOFILE (1024 on plenty of hosts) trips EMFILE on the
                # larger models. The server retries those with a plain read, so
                # this gate has to as well, otherwise it refuses to start on
                # models the server itself would have loaded quite happily.
                if e.errno != errno.EMFILE:
                    raise
                logging.warning("  too many open files, retrying without mmap: %s",
                                model.target)
                joblib.load(target, mmap_mode=None)
            logging.info("  ok: %s", model.target)
            if model.target in manifest:
                manifest[model.target]["verified_sklearn"] = sklearn.__version__
        except Exception as e:
            logging.error("  FAILED: %s (%s: %s)", model.target, type(e).__name__, e)
            broken.append(model.target)

    if broken:
        raise SystemExit(
            f"{len(broken)} model(s) could not be loaded under scikit-learn "
            f"{sklearn.__version__}: {', '.join(broken)}.\n"
            f"See the errors above for the reason. An unpickling error means "
            f"these models and this scikit-learn do not go together: either the "
            f"pinned ZENODO_RECORD is wrong for this environment, or the "
            f"scikit-learn pin moved (see app/server-requirements.txt). Refusing "
            f"to let the server start on models it cannot read."
        )


def main() -> int:
    logging.basicConfig(
        level=os.getenv("LOG_LEVEL", "INFO"),
        format="%(asctime)s %(levelname)s %(message)s",
    )

    model_dir = os.getenv("MODEL_DIR") or default_model_dir()
    record = (os.getenv("ZENODO_RECORD") or DEFAULT_ZENODO_RECORD).strip()
    wanted = selected_models(os.getenv("PARAS_MODELS", "all"))
    keep_existing = env_flag("PARAS_MODELS_KEEP_EXISTING")

    os.makedirs(model_dir, exist_ok=True)
    if not os.access(model_dir, os.W_OK):
        raise SystemExit(
            f"{model_dir} is not writable. The init service needs it mounted "
            f"read-write; only the server mounts it read-only."
        )

    logging.info("Model directory: %s", model_dir)
    logging.info("Pinned to Zenodo record %s", record)
    logging.info("Models wanted: %s", ", ".join(m.key for m in wanted))

    manifest = read_manifest(model_dir)
    missing = [m for m in wanted
               if not present_and_pinned(m, model_dir, record, manifest, keep_existing)]

    if missing:
        fetch_missing(missing, model_dir, record, manifest)
    else:
        logging.info("All %d requested model(s) present.", len(wanted))

    verify(wanted, model_dir, manifest)
    write_manifest(model_dir, manifest)
    return 0


def fetch_missing(missing: list[ModelFile], model_dir: str, record: str,
                  manifest: dict[str, dict]) -> None:
    """Download and expand the models that aren't already on disk."""
    needed = sum(m.expanded for m in missing)
    free = shutil.disk_usage(model_dir).free
    logging.info("%d model(s) to fetch: about %.2f GB to write, %.1f GB free.",
                 len(missing), needed / 1e9, free / 1e9)
    if free < needed * 1.1:
        raise SystemExit(
            f"Not enough room on {model_dir}: need about {needed / 1e9:.2f} GB, "
            f"{free / 1e9:.1f} GB free. Free up space, or set PARAS_MODELS to a "
            f"smaller set (PARAS (all substrates) alone accounts for ~3.5 GB)."
        )

    # Only talk to Zenodo once something actually needs fetching, so the common
    # case (everything already on disk) needs no network at all.
    record_files = fetch_record(record)
    for model in missing:
        manifest[model.target] = ensure(model, model_dir, record, record_files)
        write_manifest(model_dir, manifest)
    logging.info("Fetched %d model(s).", len(missing))


if __name__ == "__main__":
    sys.exit(main())
