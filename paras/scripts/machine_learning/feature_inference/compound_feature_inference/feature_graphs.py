import os
from sys import argv
from dataclasses import dataclass
from typing import Dict, List
import os


from pikachu.fingerprinting.ecfp_4 import build_ecfp_bitvector
from pikachu.general import read_smiles

import paras.data.compound_data
from paras.scripts.parsers.tabular import Tabular
from joblib import dump


FINGERPRINTS = os.path.join(os.path.dirname(paras.data.compound_data.__file__), 'fingerprints.txt')


class FeatureGraphs:
    def __init__(self):
        self.feature_graphs = []
        fingerprint_data = Tabular(FINGERPRINTS, [0])
        bits = map(int, fingerprint_data.categories[2:])
        structures = []
        compound_ids = []
        for data_id in fingerprint_data.data:
            compound_id = fingerprint_data.get_value(data_id, 'substrate_name')
            compound_smiles = fingerprint_data.get_value(data_id, 'smiles')
            structures.append(read_smiles(compound_smiles))
            compound_ids.append(compound_id)

        _, self.bit_to_substructure, _ = build_ecfp_bitvector(structures)

    def save(self, out_folder):
        if not os.path.exists(out_folder):
            os.mkdir(out_folder)

        for bit, substructure in self.bit_to_substructure.items():
            out_file = os.path.join(out_folder, f"{bit}.svg")
            substructure.draw(out_file)



if __name__ == "__main__":
    compound_features = FeatureGraphs()
    compound_features.save(argv[1])
