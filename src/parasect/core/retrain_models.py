
from joblib import dump
import os

import sklearn
from sqlalchemy import create_engine
from sqlalchemy.orm import Session

from parasect.core.constants import FINGERPRINTS_FILE, BACTERIAL_FINGERPRINTS_FILE
from parasect.model_training.rf.train_rf import train_paras_signatures, train_parasect_signatures
from parasect.model_training.train_test_splits.substrate_selection import SubstrateSelectionMode
from parasect.core.constants import INCLUDED_SUBSTRATES_FILE_PARASECT, \
    INCLUDED_SUBSTRATES_FILE_PARASECT_BACTERIAL, MODEL_METADATA_FILE, DATABASE_FILE, \
    TRAIN_DOMAINS_FILE_PARAS_ALL_SUBSTRATES, TRAIN_DOMAINS_FILE_PARASECT_BACTERIAL, TRAIN_DOMAINS_FILE_PARASECT, \
    TRAIN_DOMAINS_FILE_PARAS, INCLUDED_SUBSTRATES_FILE_PARAS, INCLUDED_SUBSTRATES_FILE_PARAS_ALL_SUBSTRATES
from parasect.core.parsing import parse_list, parse_model_metadata_file
from parasect.database.query_database import get_substrates_from_name
from parasect.core.chem import fingerprint_to_bitvector
from parasect.core.models import ModelType, Model
from parasect.core.writers import write_model_metadata_file


def model_needs_retraining(model_metadata_file: str, model_type: ModelType) -> bool:
    if not model_type.bit_count() == 1:
        raise ValueError("Must specify a single model to determine if retraining is necessary")

    model_to_version = parse_model_metadata_file(model_metadata_file)
    version = model_to_version[model_type]
    if version == sklearn.__version__:
        return False
    else:
        return True


def update_fingerprints(session: Session, hashes: list[int], included_substrates_file: str,
                        fingerprints_out_file: str) -> None:
    included_substrate_names = parse_list(included_substrates_file)
    included_substrates = []
    for name in included_substrate_names:
        included_substrates.append(get_substrates_from_name(session, name)[0])

    with open(fingerprints_out_file, 'w') as fingerprints:
        fingerprints.write("substrate_name\tsmiles")
        for fp_hash in hashes:
            fingerprints.write(f"\t{fp_hash}")
        fingerprints.write('\n')
        for substrate in included_substrates:
            fingerprints.write(f"{substrate.name}\t{substrate.smiles}")
            bitvector = fingerprint_to_bitvector(hashes, set(substrate.fingerprint))
            for value in bitvector:
                fingerprints.write(f"\t{value}")
            fingerprints.write('\n')


def update_metadata_file(model_type: ModelType, model_metadata_file: str) -> None:
    model_to_version = parse_model_metadata_file(model_metadata_file)
    model_to_version[model_type] = sklearn.__version__
    write_model_metadata_file(model_to_version, model_metadata_file)


def retrain_model(model_to_retrain: ModelType) -> Model:

    engine = create_engine(f"sqlite:///{DATABASE_FILE}")
    with Session(engine) as session:

        if model_to_retrain == ModelType.PARAS:
            model = train_paras_signatures(session, TRAIN_DOMAINS_FILE_PARAS, SubstrateSelectionMode.FIRST_ONLY,
                                           INCLUDED_SUBSTRATES_FILE_PARAS)
            return Model(model, "model.paras.gz")
        elif model_to_retrain == ModelType.PARAS_ALL_SUBSTRATES:
            model = train_paras_signatures(session, TRAIN_DOMAINS_FILE_PARAS_ALL_SUBSTRATES,
                                           SubstrateSelectionMode.FIRST_ONLY,
                                           INCLUDED_SUBSTRATES_FILE_PARAS_ALL_SUBSTRATES)
            return Model(model, "all_substrates_model.paras.gz")
        elif model_to_retrain == ModelType.PARASECT :
            model, hashes = train_parasect_signatures(session, TRAIN_DOMAINS_FILE_PARASECT,
                                                 INCLUDED_SUBSTRATES_FILE_PARASECT)
            update_fingerprints(session, hashes, INCLUDED_SUBSTRATES_FILE_PARASECT, FINGERPRINTS_FILE)
            return Model(model, "model.parasect.gz")
        elif model_to_retrain == ModelType.PARASECT_BACTERIAL:
            model, hashes = train_parasect_signatures(session, TRAIN_DOMAINS_FILE_PARASECT_BACTERIAL,
                                                 INCLUDED_SUBSTRATES_FILE_PARASECT_BACTERIAL)
            update_fingerprints(session, hashes, INCLUDED_SUBSTRATES_FILE_PARASECT_BACTERIAL,
                                BACTERIAL_FINGERPRINTS_FILE)
            return Model(model, "bacterial_model.parasect.gz")
        else:
            raise ValueError(f"Unsupported model type: {model_to_retrain}. Please attempt to retrain a single model.")


