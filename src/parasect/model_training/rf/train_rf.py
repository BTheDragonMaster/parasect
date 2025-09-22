from typing import Optional
from argparse import ArgumentParser, Namespace
from joblib import dump
from numpy.typing import NDArray
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sqlalchemy.orm import Session
from sqlalchemy import create_engine
from imblearn.under_sampling import RandomUnderSampler
import os

from parasect.core.chem import get_fingerprint_hashes, fingerprint_to_bitvector
from parasect.core.featurisation import get_domain_features
from parasect.core.writers import write_list

from parasect.database.query_database import get_domains_from_synonym, get_substrates_from_name
from parasect.core.parsing import parse_list, parse_pcs
from parasect.model_training.train_test_splits.substrate_selection import SubstrateSelectionMode, \
    map_domains_to_substrates
from parasect.model_training.train_test_splits.domain_scope import DomainScope
from parasect.model_training.rf.test_rf import test_paras_signatures, test_paras_esm, \
    test_parasect_signatures, test_parasect_esm
from parasect.database.flatfiles_from_db import write_domain_names, write_substrate_names


def parse_arguments() -> Namespace:
    """Parse arguments from command line

    :return: Arguments
    :rtype: Namespace
    """
    parser = ArgumentParser(description="Split domains into train and test set based on taxonomy")

    parser.add_argument("-db", "--database", type=str, required=True,
                        help="Path to PARASECT database")
    parser.add_argument("-o", "--out", required=True, type=str,
                        help="Output directory")
    parser.add_argument("-e", "--esm_embeddings", type=str, default=None,
                        help="Path to PCA of ESM embeddings. If given, train in ESM PCA embedding mode")
    parser.add_argument("-p", "--parasect", action="store_true",
                        help="If given, do crossvalidation for PARASECT models instead of PARAS models")
    parser.add_argument("-n", "--n_components", type=int, default=100,
                        help="Number of ESM principal components to use for domain featurisation")
    parser.add_argument('-b', "--bitvector_size", type=int, default=1024,
                        help="Size of ECFP4 bitvectors for substrate featurisation")
    parser.add_argument('-train', "--train", default=None,
                        help="List of domains to train on. If not given, train on all domains")
    parser.add_argument('-s', '--selection_mode', default="FIRST_ONLY",
                        help="Substrate selection mode for determining substrate labels")
    parser.add_argument("-i", "--included_substrates", default=None, type=str,
                        help="File containing list of included substrates. If not given, use all substrates")
    parser.add_argument('-test', "--test", default=None,
                        help="If given, test model on this domain list and write results to output folder")

    args = parser.parse_args()

    return args


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
                    included_substrates_file: str, pca_file: str, n_components: int = 100,
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
    :param n_components: number of principal components to use as features
    :type n_components: int
    :param out_path: if given, save model to out path
    :type out_path: Optional[str

    """
    domain_names = parse_list(domain_list)

    domains_pca, pca_data = parse_pcs(pca_file)
    domain_to_pcs = {
        domain: pca_data[i][:n_components]
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

    print("Training PARAS classifier...")

    classifier = train_random_forest(features, labels, out_path=out_path)

    return classifier


def train_parasect_signatures(session: Session, domain_list: str,
                              included_substrates_file: str,
                              bitvector_size: int = 1024,
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

    print("Sampling...")

    features, labels = undersampler.fit_resample(train_x, train_y)

    print("Training PARASECT classifier...")

    classifier = train_random_forest(features, labels, out_path=out_path)

    return classifier, hashes


def train_parasect_esm(session: Session, domain_list: str,
                       included_substrates_file: str,
                       pca_file: str,
                       n_components: int = 100,
                       bitvector_size: int = 256,
                       out_path: Optional[str] = None) -> tuple[RandomForestClassifier, list[int]]:
    """Train PARASECT on ESM features

    :param session: database session
    :type session: Session
    :param included_substrates_file: path to file containing included substrates
    :type included_substrates_file: str
    :param pca_file: path to pca file
    :type pca_file: str
    :param n_components: number of principal components to use for model training
    :type n_components: int
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
        domain: pca_data[i][:n_components]
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
            train_x.append(np.array(list(domain_features) + list(substrate_features)))
            if substrate in domain_to_substrates[domain]:
                train_y.append(1)
            else:
                train_y.append(0)

    train_x = np.array(train_x)
    train_y = np.array(train_y)

    undersampler = RandomUnderSampler(random_state=25051989)

    print("Sampling...")

    features, labels = undersampler.fit_resample(train_x, train_y)

    print("Training PARASECT classifier...")

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

    print("Training PARAS classifier...")

    classifier = train_random_forest(features, labels, out_path=out_path)

    return classifier


def main():
    args = parse_arguments()

    selection_mode = SubstrateSelectionMode[args.selection_mode.upper()]

    engine = create_engine(f"sqlite:///{args.database}")
    if not os.path.exists(args.out):
        os.mkdir(args.out)

    if args.parasect:
        model_path = os.path.join(args.out, "model.parasect")
    else:
        model_path = os.path.join(args.out, "model.paras")

    with Session(engine) as session:

        train_file = os.path.join(args.out, "train_domains.txt")
        if args.train is None:
            write_domain_names(session, train_file)
        else:
            write_list(parse_list(args.train), train_file)

        included_substrates_file = os.path.join(args.out, "included_substrates.txt")
        if args.included_substrates is None:
            write_substrate_names(session, included_substrates_file)
        else:
            write_list(parse_list(args.included_substrates), included_substrates_file)

        hashes = []
        if args.esm_embeddings is None:
            if not args.parasect:
                model = train_paras_signatures(session, train_file, selection_mode,
                                               included_substrates_file, out_path=model_path)
            else:
                model, hashes = train_parasect_signatures(session, train_file, included_substrates_file,
                                                          args.bitvector_size, out_path=model_path)

                print(f"Number of hashes: {len(hashes)}")
        else:
            if not args.parasect:
                model = train_paras_esm(session, train_file, selection_mode,
                                        included_substrates_file,
                                        args.esm_embeddings, args.n_components, out_path=model_path)
            else:
                model, hashes = train_parasect_esm(session, train_file, included_substrates_file, args.esm_embeddings,
                                                   args.n_components, args.bitvector_size, out_path=model_path)

                print(f"Number of hashes: {len(hashes)}")

        if hashes:
            write_list(hashes, os.path.join(args.out, "hashes.txt"), sort=False)

        if args.test:
            test_dir = os.path.join(args.out, f"test_all")
            test_dir_fungal = os.path.join(args.out, f"test_fungal")
            test_dir_bacterial = os.path.join(args.out, f"test_bacterial")

            domain_names = parse_list(args.test)
            domains = []
            write_list(domain_names, os.path.join(args.out, "test.txt"))
            for domain in domain_names:
                synonym = domain.split('|')[0]
                domains.append(get_domains_from_synonym(session, synonym)[0])

            fungal_domains = list(set(DomainScope.get_domains(session, DomainScope.FUNGAL_ONLY)).intersection(set(domains)))
            bacterial_domains = list(
                set(DomainScope.get_domains(session, DomainScope.BACTERIAL_ONLY)).intersection(set(domains)))

            print(f"\nTesting all domains...")

            if args.esm_embeddings is None:
                if not args.parasect:
                    test_paras_signatures(model, domains, included_substrates_file, test_dir)
                else:
                    assert selection_mode == SubstrateSelectionMode.ALL
                    test_parasect_signatures(session, model, domains, included_substrates_file, test_dir, hashes)

            else:
                if not args.parasect:
                    test_paras_esm(model, domains, args.esm_embeddings, included_substrates_file, test_dir,
                                   n_components=args.n_components)
                else:
                    test_parasect_esm(session, model, domains, included_substrates_file, args.esm_embeddings, test_dir,
                                      hashes, n_components=args.n_components)

            if fungal_domains and set(fungal_domains) != set(domains):
                write_list([f.get_name() for f in fungal_domains], os.path.join(args.out, "fungal_test.txt"))
                print(f"\nTesting fungal domains...")
                if args.esm_embeddings is None:
                    if not args.parasect:
                        test_paras_signatures(model, fungal_domains, included_substrates_file, test_dir_fungal)
                    else:
                        test_parasect_signatures(session, model, fungal_domains, included_substrates_file,
                                                 test_dir_fungal, hashes)
                else:
                    if not args.parasect:
                        test_paras_esm(model, fungal_domains, args.esm_embeddings, included_substrates_file,
                                       test_dir_fungal)
                    else:
                        test_parasect_esm(session, model, fungal_domains, included_substrates_file, args.esm_embeddings,
                                          test_dir_fungal, hashes, n_components=args.n_components)
            if bacterial_domains and set(bacterial_domains) != set(domains):
                write_list([b.get_name() for b in bacterial_domains], os.path.join(args.out, "bacterial_test.txt"))
                print(f"\nTesting bacterial domains...")
                if args.esm_embeddings is None:
                    if not args.parasect:
                        test_paras_signatures(model, bacterial_domains, included_substrates_file, test_dir_bacterial)
                    else:
                        test_parasect_signatures(session, model, bacterial_domains, included_substrates_file,
                                                 test_dir_bacterial, hashes)
                else:
                    if not args.parasect:
                        test_paras_esm(model, bacterial_domains, args.esm_embeddings, included_substrates_file,
                                       test_dir_bacterial)
                    else:

                        test_parasect_esm(session, model, bacterial_domains, included_substrates_file,
                                          args.esm_embeddings,
                                          test_dir_bacterial, hashes, n_components=args.n_components)


if __name__ == "__main__":
    main()
