from typing import Optional

from joblib import dump
from numpy.typing import NDArray
import numpy as np
from sklearn.ensemble import RandomForestClassifier

from parasect.core.parsing import parse_pcs, parse_fasta_file, parse_parasect_data
from parasect.core.featurisation import get_domain_features


def train_random_forest(features: NDArray[np.float64], labels: list[str], out_path: Optional[str] = None) -> None:
    """

    :param features: feature matrix
    :type features: NDArray[np.float64]
    :param labels: prediction labels
    :type labels: list[str]
    :param out_path: If given store model to specified output path
    :type out_path: Optional[str]
    """
    model = RandomForestClassifier(
        n_estimators=1000,
        n_jobs=1,
        oob_score=True,
        random_state=25051989,
        class_weight="balanced",
    )

    model.fit(features, labels, sample_weight="balanced")
    if out_path is not None:
        dump(model, out_path)


def train_rf_esm_pca():
    pass


def train_rf_signatures(signatures_fasta: str, domain_to_substrate_mapping: str, out_file: str) -> None:
    """Train RF on amino acid signatures

    :param signatures_fasta: path to fasta file containing amino acid signatures
    :type signatures_fasta: str
    :param domain_to_substrate_mapping: path to file containing domain to substrate mapping
    :type domain_to_substrate_mapping: str
    :param out_file: path to model output file
    :type out_file: str
    """
    pass

