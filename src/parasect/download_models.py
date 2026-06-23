import os
from shutil import copy
from logging import Logger

from parasect.core.models import ModelType
from parasect.core.constants import MODEL_METADATA_FILE
from parasect.core.retrain_models import retrain_model, model_needs_retraining, update_metadata_file
from parasect.core.helpers import download_and_unpack_or_fetch


def prepare_model(model_type: ModelType, model_dir: str, logger: Logger) -> str:
    """Download or retrain PARAS/PARASECT model"""
    metadata_path = os.path.join(model_dir, "model_metadata.txt")
    if not os.path.exists(metadata_path):
        copy(MODEL_METADATA_FILE, metadata_path)

    if model_needs_retraining(metadata_path, model_type):
        print("Found incompatible version of scikit-learn. Retraining..")
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
