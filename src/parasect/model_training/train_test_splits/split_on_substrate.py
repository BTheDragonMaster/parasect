from argparse import ArgumentParser, Namespace
import os

from sklearn.model_selection import train_test_split
from iterstrat.ml_stratifiers import MultilabelStratifiedShuffleSplit
from sqlalchemy.orm import Session
from sqlalchemy import create_engine
import numpy as np

from parasect.model_training.train_test_splits.domain_scope import DomainScope
from parasect.model_training.train_test_splits.substrate_selection import SubstrateSelectionMode, count_substrates, \
    count_first_valid_substrate, map_domains_to_substrates
from parasect.model_training.train_test_splits.multilabel_stratification import binarise_data
from parasect.core.writers import write_list
from parasect.model_training.train_test_splits.crossvalidation import make_crossval_sets


def parse_arguments() -> Namespace:
    """Parse arguments from command line

    :return: Arguments
    :rtype: Namespace
    """
    parser = ArgumentParser(description="Split domains into train and test set based on taxonomy")

    parser.add_argument("-db", "--database", type=str, required=True,
                        help="Path to PARASECT database")
    parser.add_argument("-o", "--output", required=True, type=str,
                        help="Output directory")
    parser.add_argument("-c", '--cutoff', default=6, type=int,
                        help="Minimum substrate count for inclusion")
    parser.add_argument("-r", '--test_ratio', default=0.25, type=float,
                        help="Target test set size")
    parser.add_argument('-f', "--fold_cross_validation", default=3, type=int,
                        help="Fold crossvalidation")

    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument('--first_only', action='store_true',
                            help="Count only the first substrate per domain")
    mode_group.add_argument('--first_valid', action='store_true',
                            help="Count only the first valid substrate per domain")

    domain_group = parser.add_mutually_exclusive_group()
    domain_group.add_argument('--fungal_only', action='store_true',
                              help="Only consider fungal domains")
    domain_group.add_argument('--bacterial_only', action='store_true',
                              help="Only consider bacterial domains")

    arguments = parser.parse_args()

    return arguments


def split_on_substrate(session: Session,
                       out_dir: str,
                       selection_mode: SubstrateSelectionMode = SubstrateSelectionMode.ALL,
                       included_domains: DomainScope = DomainScope.ALL,
                       cutoff: int = 6,
                       test_size: float = 0.25,
                       fold_cross_validation: int = 3):

    random_seed = 10012025
    domains = DomainScope.get_domains(session, included_domains)

    if selection_mode == SubstrateSelectionMode.FIRST_VALID:
        substrate_to_count = count_first_valid_substrate(session, cutoff, included_domains)
    elif selection_mode == SubstrateSelectionMode.FIRST_ONLY:
        substrate_to_count = count_substrates(session, first_only=True, cutoff=cutoff,
                                              included_domains=included_domains)
    elif selection_mode == SubstrateSelectionMode.ALL:

        substrate_to_count = count_substrates(session, first_only=False, cutoff=cutoff,
                                              included_domains=included_domains)
    else:
        raise ValueError(f"Unknown substrate selection mode: {selection_mode}")

    included_substrates = set(substrate_to_count.keys())
    domain_to_substrates = map_domains_to_substrates(domains, included_substrates=included_substrates,
                                                     selection_mode=selection_mode)
    domain_to_substrates = {
        domain: substrate
        for domain, substrate in domain_to_substrates.items()
        if substrate is not None
    }

    domains = list(domain_to_substrates.keys())

    if selection_mode == SubstrateSelectionMode.ALL:
        label_sets = []
        for domain in domains:
            substrates = domain_to_substrates[domain]
            label_sets.append(substrates)

        binary_labels, labels = binarise_data(label_sets)

        train_test_splitter = MultilabelStratifiedShuffleSplit(n_splits=1, test_size=test_size, random_state=25051989)
        train = []
        test = []
        train_y = []

        for train_x, test_x in train_test_splitter.split(domains, binary_labels):

            for idx in train_x:
                train.append(domains[idx])
                train_y.append(binary_labels[idx])

            for idx in test_x:
                test.append(domains[idx])

        make_crossval_sets(train, [domain_to_substrates[d] for d in train], out_dir, binary_labels=train_y,
                           fold_validation=fold_cross_validation, selection_mode=selection_mode)

    else:
        labels = []
        for domain in domains:
            substrate = domain_to_substrates[domain]
            labels.append(substrate)

        labels = np.array(labels)

        train, test = train_test_split(domains, stratify=labels, test_size=test_size, random_state=random_seed)
        make_crossval_sets(train, [domain_to_substrates[d] for d in train], out_dir,
                           fold_validation=fold_cross_validation, selection_mode=selection_mode)

    out_train = os.path.join(out_dir, 'train.txt')
    out_test = os.path.join(out_dir, 'test.txt')

    write_list([d.get_name() for d in train], out_train)
    write_list([d.get_name() for d in test], out_test)

    substrate_counts = {s: {"train": 0, "test": 0} for s in included_substrates}

    for domain, substrates in domain_to_substrates.items():
        if selection_mode == SubstrateSelectionMode.ALL:
            for substrate in substrates:
                if domain in train:
                    substrate_counts[substrate]["train"] += 1
                elif domain in test:
                    substrate_counts[substrate]["test"] += 1
                else:
                    raise ValueError(f"Unassigned domain: {domain.get_name()}")
        else:
            if domain in train:
                substrate_counts[substrates]["train"] += 1
            elif domain in test:
                substrate_counts[substrates]["test"] += 1
            else:
                raise ValueError(f"Unassigned domain: {domain.get_name()}")

    substrate_splits = os.path.join(out_dir, "substrate_splits.txt")

    with open(substrate_splits, 'w') as substrates_out:
        substrates_out.write("substrate_name\ttrain\ttest\n")
        for substrate, counts in substrate_counts.items():
            substrates_out.write(f"{substrate.name}\t{counts['train']}\t{counts['test']}\n")

    included_substrates_file = os.path.join(out_dir, "included_substrates.txt")
    included_substrates = list(substrate_counts.keys())
    included_substrates.sort()
    write_list(included_substrates, included_substrates_file)


def main():
    args = parse_arguments()

    if args.first_only:
        selection_mode = SubstrateSelectionMode.FIRST_ONLY
    elif args.first_valid:
        selection_mode = SubstrateSelectionMode.FIRST_VALID
    else:
        selection_mode = SubstrateSelectionMode.ALL

    if args.fungal_only:
        included_domains = DomainScope.FUNGAL_ONLY
    elif args.bacterial_only:
        included_domains = DomainScope.BACTERIAL_ONLY
    else:
        included_domains = DomainScope.ALL

    if not os.path.exists(args.output):
        os.mkdir(args.output)

    engine = create_engine(f"sqlite:///{args.database}")

    with Session(engine) as session:
        split_on_substrate(session, args.output, selection_mode=selection_mode, included_domains=included_domains,
                           test_size=args.test_ratio, cutoff=args.cutoff,
                           fold_cross_validation=args.fold_cross_validation)


if __name__ == "__main__":
    main()
