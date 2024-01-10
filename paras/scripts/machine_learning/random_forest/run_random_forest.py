import argparse
import os
from collections import OrderedDict

from paras.scripts.machine_learning.random_forest.random_forest import RandomForest, Dataset
from paras.scripts.parsers.parsers import parse_specificities, parse_morgan_fingerprint, parse_pca_file, \
    parse_pocket_features, parse_domain_list, parse_substrate_list, parse_proteinbert_features
from paras.scripts.parsers.fasta import read_fasta
from paras.scripts.feature_extraction.sequence_feature_extraction.seq_to_features import get_sequence_features_bulk
from paras.scripts.parsers.iterate_over_dir import find_crossval_pairs


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('-data', type=str, default=None, help="Path to parasect data.")
    parser.add_argument('-o', type=str, required=True, help="Path to output directory (crossval and test mode) or \
    file (train mode)")

    domain_input = parser.add_mutually_exclusive_group(required=True)
    domain_input.add_argument('-train_data', type=str,
                              help="Path to file containing list of domains. Only used in 'train' mode.")
    domain_input.add_argument('-test_data', type=str,
                              help="Path to file containing list of domains. Only used in 'test' mode.")
    domain_input.add_argument('-crossval_data', type=str,
                              help="Path to directory crossvalidation sets. Only used in 'crossval' mode.")

    parser.add_argument('-mode', type=str, default='crossval', help="Must be 'test', 'train', or 'crossval'")
    parser.add_argument('-type', type=str, required=True,
                        help="'tandem' or 'single_label'.")

    parser.add_argument('-substrates', type=str, default=None, help="Path to file containing considered substrates.")

    parser.add_argument('-f', type=str, default=None, help="Path to file containing active site signatures")
    parser.add_argument('-pca', type=str, default=None, help="Path to file containing PCs of data.")
    parser.add_argument('-pocket', type=str, default=None, help="Path to pocket file.")

    parser.add_argument('-morgan', type=str, default=None,
                        help="Path to file containing morgan fingerprints of substrates.")
    parser.add_argument('-bert', type=str, default=None, help="Path to file containing ProteinBert features.")

    parser.add_argument('-trees', type=int, default=1000, help="Nr of trees in RF")
    parser.add_argument('-threads', type=int, default=None, help="Nr of threads. -1 uses all available cores")

    parser.add_argument('-model', type=str, default=None, help="Path to classifier. Only used in 'test' mode.")

    parser.add_argument('-sampling', type=str, default=None, help="balanced, balanced_subsample, over_sample, under_sample or None")
    parser.add_argument('-one_hot', action='store_true', help="Use one-hot encoding instead of physicochemical properties for featurisation.")

    args = parser.parse_args()

    return args


def run():
    args = parse_arguments()
    if args.mode != 'train' and not os.path.exists(args.o):
        os.mkdir(args.o)

    assert args.f or args.pocket or args.pca or args.bert

    included_substrates = parse_substrate_list(args.substrates)
    domain_to_substrates = parse_specificities(args.data)

    domain_mapping_to_categories = OrderedDict()
    domain_mappings = OrderedDict()

    if args.f:
        label = "fasta"
        domain_to_seq = read_fasta(args.f)
        for domain, seq in domain_to_seq.items():
            domain_to_seq[domain] = seq.replace('X', '-')
        domain_to_features, categories = get_sequence_features_bulk(domain_to_seq, args.one_hot)
        domain_mappings[label] = domain_to_features
        domain_mapping_to_categories[label] = categories
    if args.pca:
        label = "pca"
        domain_to_pcs, categories = parse_pca_file(args.pca, return_categories=True)
        domain_mappings[label] = domain_to_pcs
        domain_mapping_to_categories[label] = categories
    if args.pocket:
        label = "pocket"
        domain_to_pocket, categories = parse_pocket_features(args.pocket, return_categories=True)
        domain_mappings[label] = domain_to_pocket
        domain_mapping_to_categories[label] = categories
    if args.bert:
        label = "bert"
        domain_to_bert, categories = parse_proteinbert_features(args.bert, return_categories=True)
        domain_mappings[label] = domain_to_bert
        domain_mapping_to_categories[label] = categories

    compound_mappings = None
    compound_mapping_to_categories = None

    if args.type == 'tandem':
        assert args.morgan

        compound_mappings = OrderedDict()
        compound_mapping_to_categories = OrderedDict()

        if args.morgan:
            label = "morgan"
            substrate_to_fingerprint, categories = parse_morgan_fingerprint(args.morgan, return_categories=True)
            compound_mappings[label] = substrate_to_fingerprint
            compound_mapping_to_categories[label] = categories

    if args.mode == 'train':
        assert args.train_data
        train_domains = parse_domain_list(args.train_data)
        rf = RandomForest(sampling_method=args.sampling, nr_trees=args.trees,
                          threads=args.threads)

        dataset = Dataset(train_domains, train_domains, [], included_substrates, domain_to_substrates,
                          domain_mappings, domain_mapping_to_categories, compound_mappings,
                          compound_mapping_to_categories, mode=args.type)
        rf.train(dataset, args.o)

    elif args.mode == 'test':
        assert args.test_data
        assert args.model
        test_domains = parse_domain_list(args.test_data)
        rf = RandomForest(model=args.model)
        dataset = Dataset(test_domains, [], test_domains, included_substrates, domain_to_substrates,
                          domain_mappings, domain_mapping_to_categories, compound_mappings,
                          compound_mapping_to_categories, mode=args.type)
        rf.record_test_results(dataset, args.o, mode='test')

    elif args.mode == 'crossval':
        assert args.crossval_data
        crossval_pairs = find_crossval_pairs(args.crossval_data)
        for i, crossval_pair in enumerate(crossval_pairs):
            crossval_out = os.path.join(args.o, f"crossval_{i + 1}")
            if not os.path.exists(crossval_out):
                os.mkdir(crossval_out)

            train_domain_path, test_domain_path = crossval_pair
            train_domains = parse_domain_list(train_domain_path)
            test_domains = parse_domain_list(test_domain_path)
            all_domains = train_domains + test_domains

            rf = RandomForest(sampling_method=args.sampling, nr_trees=args.trees,
                              threads=args.threads)

            dataset = Dataset(all_domains, train_domains, test_domains, included_substrates, domain_to_substrates,
                              domain_mappings, domain_mapping_to_categories, compound_mappings,
                              compound_mapping_to_categories, mode=args.type)

            rf.train(dataset)
            rf.record_test_results(dataset, crossval_out, mode='both')
            rf.write_feature_importances(dataset.categories, crossval_out)


if __name__ == "__main__":
    run()
