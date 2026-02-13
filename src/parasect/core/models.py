from dataclasses import dataclass
import os
from joblib import dump
from enum import IntFlag

from sklearn.ensemble import RandomForestClassifier


class ModelType(IntFlag):
    PARAS = 1
    PARASECT = 2
    PARAS_ALL_SUBSTRATES = 4
    PARASECT_BACTERIAL = 8
    ANTISMASH_MODELS = PARAS_ALL_SUBSTRATES | PARASECT_BACTERIAL
    ALL_MODELS = PARAS | PARASECT | PARAS_ALL_SUBSTRATES | PARASECT_BACTERIAL


@dataclass
class Model:
    model: RandomForestClassifier
    file_name: str

    def save(self, out_dir: str) -> None:
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        dump(self.model, os.path.join(out_dir, self.file_name))
