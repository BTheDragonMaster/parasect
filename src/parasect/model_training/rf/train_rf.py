from typing import Optional

from joblib import dump
from numpy.typing import NDArray
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sqlalchemy.orm import Session
from imblearn.under_sampling import RandomUnderSampler

from parasect.core.chem import get_fingerprint_hashes, fingerprint_to_bitvector
from parasect.core.featurisation import get_domain_features
from parasect.database.query_database import get_domains_from_synonym, get_substrates_from_name
from parasect.core.parsing import parse_list, parse_pcs
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

    model.fit(features, labels)
    if out_path is not None:
        dump(model, out_path)

    return model


def train_paras_esm(session: Session, domain_list: str, selection_mode: SubstrateSelectionMode,
                    included_substrates_file: str, pca_file: str,
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
    :param pca_file: path to file containing pca of esm embeddings
    :type pca_file: str
    :param out_path: if given, save model to out path
    :type out_path: Optional[str

    """
    domain_names = parse_list(domain_list)

    domains_pca, pca_data = parse_pcs(pca_file)
    domain_to_pcs = {
        domain: pca_data[i]
        for i, domain in enumerate(domains_pca)
        if domain in domain_names
    }

    domains = []
    for domain in domain_names:
        synonym = domain.split('|')[0]
        domains.append(get_domains_from_synonym(session, synonym)[0])

    domains.sort(key=lambda x: x.get_name())
    included_substrate_names = parse_list(included_substrates_file)
    included_substrates = set()
    for name in included_substrate_names:
        included_substrates.add(get_substrates_from_name(session, name)[0])

    domain_to_substrate = map_domains_to_substrates(domains, included_substrates, selection_mode=selection_mode)

    domain_to_substrate = {
        domain: substrate
        for domain, substrate in domain_to_substrate.items()
        if substrate is not None
    }

    features = np.array([domain_to_pcs[domain.get_name()] for domain in domains])
    labels = np.array([domain_to_substrate[domain].name for domain in domains])

    classifier = train_random_forest(features, labels, out_path=out_path)

    return classifier


def train_parasect_signatures(session: Session, domain_list: str,
                              included_substrates_file: str,
                              bitvector_size: int = 256,
                              out_path: Optional[str] = None) -> tuple[RandomForestClassifier, list[int]]:
    """Train PARASECT on extended signatures

    :param session: database session
    :type session: Session
    :param included_substrates_file: path to file containing included substrates
    :type included_substrates_file: str
    :param bitvector_size: size of bitvector
    :type bitvector_size: int
    :param domain_list: file path containing list of domains to train on
    :type domain_list: str
    :param out_path: if given, save model to out path
    :type out_path: Optional[str
    :returns: random forest classifier and sorted hashes used for bitvector construction
    :rtype: tuple[RandomForestClassifier, list[int]]

    """

    domain_names = parse_list(domain_list)

    domains = []
    for domain in domain_names:
        synonym = domain.split('|')[0]
        domains.append(get_domains_from_synonym(session, synonym)[0])

    domains.sort(key=lambda x: x.get_name())
    included_substrate_names = parse_list(included_substrates_file)
    included_substrates = set()
    for name in included_substrate_names:
        included_substrates.add(get_substrates_from_name(session, name)[0])

    hashes = get_fingerprint_hashes(included_substrates, bitvector_size)

    domain_to_substrates = map_domains_to_substrates(domains, included_substrates,
                                                     selection_mode=SubstrateSelectionMode.ALL)

    domain_to_substrates = {
        domain: substrates
        for domain, substrates in domain_to_substrates.items()
        if substrates is not None
    }
    train_x = []
    train_y = []
    for domain in domain_to_substrates:
        for substrate in included_substrates:
            domain_features = get_domain_features(domain.extended_signature)
            substrate_features = fingerprint_to_bitvector(hashes, set(substrate.fingerprint))
            train_x.append(np.array(domain_features + substrate_features))
            if substrate in domain_to_substrates[domain]:
                train_y.append(1)
            else:
                train_y.append(0)

    undersampler = RandomUnderSampler(random_state=25051989)

    features, labels = undersampler.fit_resample(train_x, train_y)

    classifier = train_random_forest(features, labels, out_path=out_path)

    return classifier, hashes


def train_parasect_esm(session: Session, domain_list: str,
                       included_substrates_file: str,
                       pca_file: str,
                       bitvector_size: int = 256,
                       out_path: Optional[str] = None) -> tuple[RandomForestClassifier, list[int]]:
    """Train PARASECT on ESM features

    :param session: database session
    :type session: Session
    :param included_substrates_file: path to file containing included substrates
    :type included_substrates_file: str
    :param pca_file: path to pca file
    :type pca_file: str
    :param bitvector_size: size of bitvector
    :type bitvector_size: int
    :param domain_list: file path containing list of domains to train on
    :type domain_list: str
    :param out_path: if given, save model to out path
    :type out_path: Optional[str
    :returns: random forest classifier and sorted hashes used for bitvector construction
    :rtype: tuple[RandomForestClassifier, list[int]]

    """

    domain_names = parse_list(domain_list)

    domains_pca, pca_data = parse_pcs(pca_file)
    domain_to_pcs = {
        domain: pca_data[i]
        for i, domain in enumerate(domains_pca)
        if domain in domain_names
    }

    domains = []
    for domain in domain_names:
        synonym = domain.split('|')[0]
        domains.append(get_domains_from_synonym(session, synonym)[0])

    domains.sort(key=lambda x: x.get_name())
    included_substrate_names = parse_list(included_substrates_file)
    included_substrates = set()
    for name in included_substrate_names:
        included_substrates.add(get_substrates_from_name(session, name)[0])

    hashes = get_fingerprint_hashes(included_substrates, bitvector_size)

    domain_to_substrates = map_domains_to_substrates(domains, included_substrates,
                                                     selection_mode=SubstrateSelectionMode.ALL)

    domain_to_substrates = {
        domain: substrates
        for domain, substrates in domain_to_substrates.items()
        if substrates is not None
    }

    train_x = []
    train_y = []
    for domain in domain_to_substrates:
        for substrate in included_substrates:
            domain_features = domain_to_pcs[domain.get_name()]
            substrate_features = fingerprint_to_bitvector(hashes, set(substrate.fingerprint))
            train_x.append(np.array(domain_features + substrate_features))
            if substrate in domain_to_substrates[domain]:
                train_y.append(1)
            else:
                train_y.append(0)

    train_x = np.array(train_x)
    train_y = np.array(train_y)

    undersampler = RandomUnderSampler(random_state=25051989)

    features, labels = undersampler.fit_resample(train_x, train_y)

    classifier = train_random_forest(features, labels, out_path=out_path)

    return classifier, hashes


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

    domains = []
    for domain in domain_names:
        synonym = domain.split('|')[0]
        domains.append(get_domains_from_synonym(session, synonym)[0])

    domains.sort(key=lambda x: x.get_name())
    included_substrate_names = parse_list(included_substrates_file)
    included_substrates = set()
    for name in included_substrate_names:
        included_substrates.add(get_substrates_from_name(session, name)[0])

    domain_to_substrate = map_domains_to_substrates(domains, included_substrates, selection_mode=selection_mode)

    domain_to_substrate = {
        domain: substrate
        for domain, substrate in domain_to_substrate.items()
        if substrate is not None
    }

    features = np.array([get_domain_features(domain.extended_signature) for domain in domains])
    labels = np.array([domain_to_substrate[domain].name for domain in domains])

    classifier = train_random_forest(features, labels, out_path=out_path)

    return classifier
