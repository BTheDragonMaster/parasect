import os
from enum import Enum
from typing import Optional

from sqlalchemy.orm import Session
from sqlalchemy import create_engine
from argparse import ArgumentParser, Namespace


from parasect.core.parsing import iterate_over_dir, parse_list
from parasect.database.query_database import get_domains_from_synonym
from parasect.model_training.rf.train_rf import train_paras_signatures, train_paras_esm
from parasect.model_training.rf.test_rf import test_paras_signatures, test_paras_esm
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
    parser.add_argument("-p", "--esm_pca", type=str, default=None,
                        help="Path to PCA of ESM embeddings. If given, train in ESM PCA embedding mode")

    args = parser.parse_args()

    return args


class TrainingMode(Enum):
    SIGNATURES = 1
    ESM_PCA = 2


def do_crossvalidation(session: Session, crossvalidation_folder: str, out_dir: str,
                       esm_embeddings: Optional[str] = None) -> None:
    """

    :param session:
    :type session:
    :param crossvalidation_folder:
    :type crossvalidation_folder:
    :param out_dir:
    :type out_dir:
    :param esm_embeddings:
    :type esm_embeddings:
    :return:
    :rtype:
    """
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    included_substrates_file = os.path.join(crossvalidation_folder, "included_substrates.txt")
    crossvalidation_train_files = []
    crossvalidation_test_files = []
    if "_valid" in crossvalidation_folder:
        selection_mode = SubstrateSelectionMode.FIRST_VALID
    elif "_first" in crossvalidation_folder:
        selection_mode = SubstrateSelectionMode.FIRST_ONLY
    else:
        selection_mode = SubstrateSelectionMode.FIRST_ONLY

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

        if esm_embeddings is None:
            model = train_paras_signatures(session, train_file, selection_mode, included_substrates_file)
        else:
            model = train_paras_esm(session, train_file, selection_mode, included_substrates_file,
                                    esm_embeddings)

        crossval_dir = os.path.join(out_dir, f"crossvalidation_{i + 1}_all")
        crossval_dir_fungal = os.path.join(out_dir, f"crossvalidation_{i + 1}_fungal")
        crossval_dir_bacterial = os.path.join(out_dir, f"crossvalidation_{i + 1}_bacterial")
        print(f"\nTesting all domains... (crossvalidation set {i + 1})")
        if esm_embeddings is None:
            test_paras_signatures(model, domains, included_substrates_file, crossval_dir)
        else:
            test_paras_esm(model, domains, esm_embeddings, included_substrates_file, crossval_dir)
        if fungal_domains:
            print(f"\nTesting fungal domains... (crossvalidation set {i + 1})")
            if esm_embeddings is None:
                test_paras_signatures(model, fungal_domains, included_substrates_file, crossval_dir_fungal)
            else:
                test_paras_esm(model, fungal_domains, esm_embeddings, included_substrates_file, crossval_dir_fungal)
        if bacterial_domains:
            print(f"\nTesting bacterial domains... (crossvalidation set {i + 1})")
            if esm_embeddings is None:
                test_paras_signatures(model, bacterial_domains, included_substrates_file, crossval_dir_bacterial)
            else:
                test_paras_esm(model, bacterial_domains, esm_embeddings, included_substrates_file, crossval_dir_bacterial)


if __name__ == "__main__":
    args = parse_arguments()

    engine = create_engine(f"sqlite:///{args.database}")
    with Session(engine) as session:
        do_crossvalidation(session, args.crossval_folder, args.output, args.esm_pca)