import tarfile
import logging
from urllib.parse import urlsplit
from urllib.request import urlopen, Request
from pathlib import Path
from typing import Optional
from parasect.core.parsing import parse_smiles_mapping
import os
from shutil import copy

from parasect.core.models import ModelType
from parasect.core.constants import MODEL_METADATA_FILE
from parasect.core.retrain_models import retrain_model, model_needs_retraining, update_metadata_file

logger = logging.getLogger(__name__)

def prepare_model(model_type: ModelType, model_dir: str) -> str:
    """Download or retrain PARAS/PARASECT model"""
    metadata_path = os.path.join(model_dir, "model_metadata.txt")
    if not os.path.exists(metadata_path):
        copy(MODEL_METADATA_FILE, metadata_path)

    if model_needs_retraining(metadata_path, model_type):
        logger.info("Found incompatible version of scikit-learn. Retraining..")
        model = retrain_model(model_type)
        model_path = os.path.join(model_dir, model.file_name)
        model.save(model_dir)
        update_metadata_file(model_type, metadata_path)

    else:
        if model_type == ModelType.PARAS_ALL_SUBSTRATES:
            model_path = download_and_unpack_or_fetch(
                r"https://zenodo.org/records/17224548/files/all_substrates_model.paras.gz?download=1",
                model_dir, logger)
        elif model_type == ModelType.PARAS:
            model_path = download_and_unpack_or_fetch(
                r"https://zenodo.org/records/17224548/files/model.paras.gz?download=1",
                model_dir, logger)
        elif model_type == ModelType.PARASECT:
            model_path = download_and_unpack_or_fetch(
                r"https://zenodo.org/records/17224548/files/model.parasect.gz?download=1",
                model_dir, logger)
        elif model_type == ModelType.PARASECT_BACTERIAL:
            model_path = download_and_unpack_or_fetch(
                r"https://zenodo.org/records/17224548/files/bacterial_model.parasect.gz?download=1",
                model_dir, logger)
        else:
            raise ValueError("Unknown model type")

    return model_path


def prepare_substrates(smiles_mapping: Optional[str]) -> tuple[Optional[list[str]], Optional[list[str]]]:
    """Return substrate names and substrate SMILES from SMILES mapping

    :param smiles_mapping: path to file containing substrate names in column 1 and SMILES strings in column 2

    :returns: list of substrate names and list substrate SMILES if SMILES mapping exists, tuple of (None,None) otherwise

    """
    if smiles_mapping is not None:
        substrates = parse_smiles_mapping(smiles_mapping)
        substrate_names = [s.name for s in substrates]
        substrate_smiles = [s.smiles for s in substrates]
    else:
        substrate_names = None
        substrate_smiles = None

    return substrate_names, substrate_smiles

def prepare_folders(out_dir: str, temp_dir: Optional[str], model_dir: Optional[str]) -> tuple[str, str]:
    """Prepare folders for output

    :param out_dir: Output directory
    :param temp_dir: Temporary directory
    :param model_dir: Model directory

    :returns: paths to temporary directory and model directory
    """

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    if temp_dir is None:
        temp_dir = os.path.join(out_dir, "temp")

    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)

    if model_dir is None:
        model_dir = os.path.join(out_dir, "model")

    if not os.path.exists(model_dir):
        os.mkdir(model_dir)

    return temp_dir, model_dir


def download_and_unpack_or_fetch(
    url: str,
    dest_dir: str,
    logger: logging.Logger
) -> str:
    """
    Download from `url` into `dest_dir` (creating it if needed).

    - Archives (.tar.gz/.tgz):
        • If already downloaded & unpacked, returns existing unpacked path.
        • Otherwise downloads & extracts, then returns unpacked folder path.
    - Plain files:
        • If already downloaded, returns existing file path.
        • Otherwise downloads it, then returns the file path.

    Logs each step and shows download progress via the given logger.

    :param url: HTTP(S) URL of the file or archive to download.
    :param dest_dir: Directory to save the downloaded file or unpacked archive.
    :param logger: Logger for status/progress messages.
    :return: Absolute path to the downloaded file or unpacked archive.
    """
    dest = Path(dest_dir).expanduser().absolute()
    logger.info(f"Using destination directory: {dest}")
    dest.mkdir(parents=True, exist_ok=True)

    filename = Path(urlsplit(url).path).name
    local_path = dest / filename
    logger.info(f"Local path will be: {local_path}")

    def _download():
        req = Request(url, headers={'User-Agent': 'Python urllib'})
        with urlopen(req) as resp, open(local_path, 'wb') as out:
            total_size = int(resp.getheader('Content-Length') or 0)
            if total_size:
                logger.info(f"Starting download of {filename} ({total_size} bytes)")
                next_pct = 10
            else:
                logger.info(f"Starting download of {filename} (size unknown)")
                threshold = 1 * 1024 * 1024  # 1 MiB
                next_bytes = threshold

            downloaded = 0
            chunk_size = 8192

            while True:
                chunk = resp.read(chunk_size)
                if not chunk:
                    break
                out.write(chunk)
                downloaded += len(chunk)

                if total_size:
                    pct = downloaded * 100 / total_size
                    if pct >= next_pct:
                        logger.info(f"Downloaded {downloaded}/{total_size} bytes ({pct:.0f}%)")
                        next_pct += 10
                else:
                    if downloaded >= next_bytes:
                        logger.info(f"Downloaded {downloaded} bytes")
                        next_bytes += threshold

            logger.info("Download complete")

    is_tar = filename.lower().endswith(('.tar.gz', '.tgz'))

    if is_tar:
        logger.info(f"Detected archive format for {filename}")

        def _get_roots(tar_path):
            with tarfile.open(tar_path, 'r:gz') as tar:
                return {Path(m.name).parts[0] for m in tar.getmembers() if m.name.strip()}

        # If archive exists, see if its contents are already unpacked
        if local_path.exists():
            logger.info(f"Archive already exists at {local_path}, checking extraction")
            try:
                roots = _get_roots(local_path)
                if len(roots) == 1:
                    candidate = dest / roots.pop()
                    if candidate.exists():
                        logger.info(f"Found existing extracted folder {candidate}, skipping.")
                        return str(candidate)
                else:
                    others = [p for p in dest.iterdir() if p != local_path]
                    if others:
                        logger.info(f"Found already-extracted contents in {dest}, skipping.")
                        return str(dest)
            except tarfile.TarError:
                logger.warning(f"Archive at {local_path} seems corrupt. Re-downloading.")

        # Download & extract
        _download()
        logger.info(f"Extracting {local_path} to {dest}")
        with tarfile.open(local_path, 'r:gz') as tar:
            roots = {Path(m.name).parts[0] for m in tar.getmembers() if m.name.strip()}
            tar.extractall(path=dest)

        if len(roots) == 1:
            unpacked = dest / roots.pop()
            logger.info(f"Extraction complete: single folder {unpacked}")
            return str(unpacked)
        else:
            logger.info(f"Extraction complete: multiple items, returning {dest}")
            return str(dest)

    else:
        # Plain file
        logger.info(f"Detected single file for {filename}")
        if not local_path.exists():
            logger.info(f"{local_path} not found locally; downloading.")
            _download()
        else:
            logger.info(f"File already exists at {local_path}, skipping download.")
        return str(local_path)
