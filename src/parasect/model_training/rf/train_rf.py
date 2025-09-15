from typing import Optional

from joblib import dump
from numpy.typing import NDArray
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sqlalchemy.orm import Session

from parasect.core.featurisation import get_domain_features
from parasect.database.query_database import get_domains_from_names, get_substrates_from_name
from parasect.core.parsing import parse_list
from parasect.model_training.train_test_splits.substrate_selection import SubstrateSelectionMode, \
    map_domains_to_substrates


def train_random_forest(features: NDArray[np.float64], labels: NDArray[np.str],
                        out_path: Optional[str] = None) -> RandomForestClassifier:
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

    return model


def train_rf_esm_pca():
    pass


def train_paras_signatures(session: Session, domain_list: str,
                           selection_mode: SubstrateSelectionMode,
                           included_substrates_file: str,
                           out_path: Optional[str] = None) -> RandomForestClassifier:
    """

    :param session: database session
    :type session: Session
    :param selection_mode: determines which substrates are considered
    :type selection_mode: SubstrateSelectionMode
    :param included_substrates_file: path to file containing included substrates
    :type included_substrates_file: str
    :param domain_list: file path containing list of domains to train on
    :type domain_list: str
    :param out_path: if given, save model to out path
    :type out_path: Optional[str

    """

    if selection_mode not in [SubstrateSelectionMode.FIRST_VALID, SubstrateSelectionMode.FIRST_ONLY]:
        raise ValueError(f"Unsupported substrate selection mode for training PARAS: {selection_mode.name}")

    domain_names = parse_list(domain_list)

    domains = get_domains_from_names(session, domain_names)

    domains.sort(key=lambda x: x.get_name())
    included_substrate_names = parse_list(included_substrates_file)
    included_substrates = set()
    for name in included_substrate_names:
        included_substrates.add(get_substrates_from_name(session, name)[0])

    domain_to_substrate = map_domains_to_substrates(domains, included_substrates, selection_mode=selection_mode)

    features = np.array([get_domain_features(domain.extended_signature) for domain in domains])
    labels = np.array([domain_to_substrate[domain].name for domain in domains])

    classifier = train_random_forest(features, labels, out_path=out_path)

    return classifier









