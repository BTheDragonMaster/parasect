import os
from enum import Enum
from typing import Optional

from sqlalchemy.orm import Session
from sqlalchemy import create_engine
from argparse import ArgumentParser, Namespace


from parasect.core.parsing import iterate_over_dir, parse_list
from parasect.database.query_database import get_domains_from_synonym
from parasect.model_training.rf.train_rf import train_paras_signatures, train_paras_esm, \
    train_parasect_esm, train_parasect_signatures
from parasect.model_training.rf.test_rf import test_paras_signatures, test_paras_esm, test_parasect_signatures, \
    test_parasect_esm
from parasect.model_training.train_test_splits.domain_scope import DomainScope
from parasect.model_training.train_test_splits.substrate_selection import SubstrateSelectionMode


def parse_arguments() -> Namespace:
    """Parse arguments from command line

    :return: Arguments
    :rtype: Namespace
    """
    parser = ArgumentParser(description="Split domains into train and test set based on taxonomy")

    parser.add_argument("-db", "--database", type=str, required=True,
                        help="Path to PARASECT database")
    parser.add_argument("-c", "--crossval_folder", type=str, required=True,
                        help="Path to crossvalidation folder")
    parser.add_argument("-o", "--output", required=True, type=str,
                        help="Output directory")
    parser.add_argument("-e", "--esm_pca", type=str, default=None,
                        help="Path to PCA of ESM embeddings. If given, train in ESM PCA embedding mode")
    parser.add_argument("-p", "--parasect", action="store_true",
                        help="If given, do crossvalidation for PARASECT models instead of PARAS models")
    parser.add_argument("-n", "--n_components", type=int, default=100,
                        help="Number of ESM principal components to use for domain featurisation")
    parser.add_argument('-b', "--bitvector_size", type=int, default=1024,
                        help="Size of ECFP4 bitvectors for substrate featurisation")

    args = parser.parse_args()

    return args


class TrainingMode(Enum):
    SIGNATURES = 1
    ESM_PCA = 2


def do_crossvalidation(session: Session, crossvalidation_folder: str, out_dir: str,
                       esm_embeddings: Optional[str] = None, parasect: bool = False,
                       bitvector_size: int = 256, n_components: int = 100) -> None:
    """

    :param session:
    :type session:
    :param crossvalidation_folder:
    :type crossvalidation_folder:
    :param out_dir:
    :type out_dir:
    :param esm_embeddings:
    :type esm_embeddings:
    :param parasect: if True, do crossvalidation for PARASECT models instead of PARAS models
    :type parasect: bool
    :param bitvector_size: Size of ECFP4 bitvectors for substrate featurisation
    :type bitvector_size: int
    :param n_components: Number of ESM principal components to use for domain featurisation
    :type n_components: int
    :return:
    :rtype:
    """
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    included_substrates_file = os.path.join(crossvalidation_folder, "parasect.txt")
    crossvalidation_train_files = []
    crossvalidation_test_files = []
    if "_valid" in crossvalidation_folder:
        selection_mode = SubstrateSelectionMode.FIRST_VALID
    elif "_first" in crossvalidation_folder:
        selection_mode = SubstrateSelectionMode.FIRST_ONLY
    else:
        selection_mode = SubstrateSelectionMode.ALL

    for file_name, file_path in iterate_over_dir(crossvalidation_folder, '.txt'):
        if file_name.startswith("test_crossvalidation"):
            crossvalidation_test_files.append(file_path)
        if file_name.startswith("train_crossvalidation"):
            crossvalidation_train_files.append(file_path)

    crossvalidation_train_files.sort()
    crossvalidation_test_files.sort()

    assert len(crossvalidation_train_files) == len(crossvalidation_test_files)

    for i, train_file in enumerate(crossvalidation_train_files):

        test_file = crossvalidation_test_files[i]

        domain_names = parse_list(test_file)
        domains = []
        for domain in domain_names:
            synonym = domain.split('|')[0]
            domains.append(get_domains_from_synonym(session, synonym)[0])

        fungal_domains = list(set(DomainScope.get_domains(session, DomainScope.FUNGAL_ONLY)).intersection(set(domains)))
        bacterial_domains = list(set(DomainScope.get_domains(session, DomainScope.BACTERIAL_ONLY)).intersection(set(domains)))

        hashes = []
        if esm_embeddings is None:
            if not parasect:
                model = train_paras_signatures(session, train_file, selection_mode, included_substrates_file)
            else:
                model, hashes = train_parasect_signatures(session, train_file, included_substrates_file,
                                                          bitvector_size)

                print(f"Number of hashes: {len(hashes)}")
        else:
            if not parasect:
                model = train_paras_esm(session, train_file, selection_mode, included_substrates_file,
                                        esm_embeddings, n_components)
            else:
                model, hashes = train_parasect_esm(session, train_file, included_substrates_file, esm_embeddings,
                                                   n_components, bitvector_size)

                print(f"Number of hashes: {len(hashes)}")

        crossval_dir = os.path.join(out_dir, f"crossvalidation_{i + 1}_all")
        crossval_dir_fungal = os.path.join(out_dir, f"crossvalidation_{i + 1}_fungal")
        crossval_dir_bacterial = os.path.join(out_dir, f"crossvalidation_{i + 1}_bacterial")

        print(f"\nTesting all domains... (crossvalidation set {i + 1})")

        if esm_embeddings is None:
            if not parasect:
                test_paras_signatures(model, domains, included_substrates_file, crossval_dir)
            else:
                assert selection_mode == SubstrateSelectionMode.ALL
                test_parasect_signatures(session, model, domains, included_substrates_file, crossval_dir, hashes)

        else:
            if not parasect:
                test_paras_esm(model, domains, esm_embeddings, included_substrates_file, crossval_dir,
                               n_components=n_components)
            else:
                test_parasect_esm(session, model, domains, included_substrates_file, esm_embeddings, crossval_dir,
                                  hashes, n_components=n_components)

        if fungal_domains and set(fungal_domains) != set(domains):
            print(f"\nTesting fungal domains... (crossvalidation set {i + 1})")
            if esm_embeddings is None:
                if not parasect:
                    test_paras_signatures(model, fungal_domains, included_substrates_file, crossval_dir_fungal)
                else:
                    test_parasect_signatures(session, model, fungal_domains, included_substrates_file,
                                             crossval_dir_fungal, hashes)
            else:
                if not parasect:
                    test_paras_esm(model, fungal_domains, esm_embeddings, included_substrates_file, crossval_dir_fungal)
                else:
                    test_parasect_esm(session, model, fungal_domains, included_substrates_file, esm_embeddings,
                                      crossval_dir_fungal, hashes, n_components=n_components)
        if bacterial_domains and set(bacterial_domains) != set(domains):
            print(f"\nTesting bacterial domains... (crossvalidation set {i + 1})")
            if esm_embeddings is None:
                if not parasect:
                    test_paras_signatures(model, bacterial_domains, included_substrates_file, crossval_dir_bacterial)
                else:
                    test_parasect_signatures(session, model, bacterial_domains, included_substrates_file,
                                             crossval_dir_bacterial, hashes)
            else:
                if not parasect:
                    test_paras_esm(model, bacterial_domains, esm_embeddings, included_substrates_file,
                                   crossval_dir_bacterial)
                else:

                    test_parasect_esm(session, model, bacterial_domains, included_substrates_file, esm_embeddings,
                                      crossval_dir_bacterial, hashes, n_components=n_components)


if __name__ == "__main__":
    args = parse_arguments()

    engine = create_engine(f"sqlite:///{args.database}")
    with Session(engine) as session:
        do_crossvalidation(session, args.crossval_folder, args.output, args.esm_pca,
                           args.parasect, args.bitvector_size, args.n_components)