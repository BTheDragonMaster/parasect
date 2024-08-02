import os
from sys import argv
import argparse

from paras.scripts.parsers.external_parsers import parse_sandpuma_data, parse_nrpspredictor_data, \
    parse_adenpredictor_data
from paras.scripts.parsers.parsers import parse_domain_list, parse_test_results, parse_specificities

import paras.data.structure_data
import paras.data.train_test_splits
import paras.data.benchmarking
import paras.data
from paras.scripts.parsers.iterate_over_dir import iterate_over_dir

benchmarking_dir = os.path.dirname(paras.data.benchmarking.__file__)


SUBSTRATE_FILE = os.path.join(os.path.dirname(paras.data.structure_data.__file__), 'included_substrates.txt')
DOMAINS_CLASS_SANDPUMA = os.path.join(benchmarking_dir, 'benchmark_class_filtered.txt')
DOMAINS_PHYLOGENY_SANDPUMA = os.path.join(benchmarking_dir, 'benchmark_phylogeny_filtered.txt')
DOMAINS_CLASS = os.path.join(os.path.dirname(paras.data.train_test_splits.__file__), 'class/test_filtered.txt')
DOMAINS_PHYLOGENY = os.path.join(os.path.dirname(paras.data.train_test_splits.__file__), 'phylogeny/test_filtered.txt')
PARAS_RESULTS_DIR = os.path.join(benchmarking_dir, 'paras_results')
PARASECT_RESULTS_DIR = os.path.join(benchmarking_dir, 'parasect_results')
PARASECT_DATASET = os.path.join(os.path.dirname(paras.data.__file__), "parasect_dataset.txt")


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', type=str, help="Path to PARAS & PARASECT results folder.")
    parser.add_argument('-o', type=str, help="Path to output file.")
    parser.add_argument('-d', type=str, help="Path to domain list representing holdout test set.")
    parser.add_argument('-m', type=str, default='class', help="Method used to create holdout test set. Must be 'class' or 'phylogeny'.")
    parser.add_argument('-s', action='store_true', help="Add sandpuma results to summary file.")
    parser.add_argument('-a', action='store_true', help="Add adenpredictor results to summary file.")
    parser.add_argument('-n', action='store_true', help="Add nrpspredictor results to summary file.")

    args = parser.parse_args()

    return args

def update_data(test_results, method, domain_to_method_to_prediction):

    domain_to_prediction = parse_test_results(test_results)

    for domain in domain_to_method_to_prediction:

        domain_to_method_to_prediction[domain][method] = domain_to_prediction[domain]


def get_summary(results_dir, out_file, holdout_set, split_method, sandpuma=True, nrpspredictor=True, adenpredictor=True):

    assert split_method in ["phylogeny", 'class']

    domain_to_spec = parse_specificities(PARASECT_DATASET)

    domain_list = parse_domain_list(holdout_set)

    domain_to_method_to_prediction = {}
    for domain in domain_list:
        domain_to_method_to_prediction[domain] = {}
        domain_to_method_to_prediction[domain]["substrate"] = '|'.join(domain_to_spec[domain])

    for model_name, model_dir in iterate_over_dir(results_dir, get_dirs=True):
        for extraction_method, extraction_dir in iterate_over_dir(model_dir, get_dirs=True):
            if extraction_method == 'hmm':
                for featurisation_name, featurisation_dir in iterate_over_dir(extraction_dir, get_dirs=True):
                    if split_method in featurisation_name:
                        featurisation_type = featurisation_name.split(f'{split_method}_')[1]
                        assert featurisation_type in ['sequence', 'structure', 'sequence_structure']
                        if featurisation_type == 'sequence_structure':
                            featurisation = "Sequence + Structure"
                        else:
                            featurisation = featurisation_type.capitalize()
                        for test_dir_name, test_dir in iterate_over_dir(featurisation_dir, get_dirs=True):
                            if test_dir_name == 'test_performance':
                                test_results = os.path.join(test_dir, 'test_results.txt')
                                update_data(test_results, f"{model_name.upper()} ({featurisation})",
                                            domain_to_method_to_prediction)

    if sandpuma:
        add_sandpuma_data(domain_to_method_to_prediction)
    if nrpspredictor:
        add_nrpspredictor_data(domain_to_method_to_prediction)
    if adenpredictor:
        add_adenpredictor_data(domain_to_method_to_prediction)

    with open(out_file, 'w') as out:
        paras_methods = []
        parasect_methods = []
        substrate_methods = []
        other_methods = []
        for method in domain_to_method_to_prediction[list(domain_to_method_to_prediction.keys())[0]]:
            if 'PARASECT' in method:
                parasect_methods.append(method)
            elif 'PARAS' in method:
                paras_methods.append(method)
            elif method == 'substrate':
                substrate_methods.append(method)
            else:
                other_methods.append(method)

        paras_methods.sort(key = lambda x: len(x))
        parasect_methods.sort(key=lambda x: len(x))
        other_methods.sort()

        methods = substrate_methods + paras_methods + parasect_methods + other_methods
        out.write('domain\t')
        out.write('\t'.join(methods))
        out.write('\n')

        for domain, method_to_prediction in domain_to_method_to_prediction.items():
            out.write(domain)
            for method in methods:
                prediction = method_to_prediction[method]
                out.write(f"\t{prediction}")
            out.write('\n')


def add_sandpuma_data(domain_to_method_to_prediction):
    domain_to_prediction = parse_sandpuma_data()

    for domain in domain_to_method_to_prediction:
        domain_to_method_to_prediction[domain]["SANDPUMA"] = '|'.join(domain_to_prediction[domain].sandpuma)
        domain_to_method_to_prediction[domain]["SANDPUMA (ensemble)"] = '|'.join(domain_to_prediction[domain].ensemble)
        domain_to_method_to_prediction[domain]["SANDPUMA (ASM)"] = '|'.join(domain_to_prediction[domain].asm)
        domain_to_method_to_prediction[domain]["SANDPUMA (SVM)"] = '|'.join(domain_to_prediction[domain].svm)
        domain_to_method_to_prediction[domain]["SANDPUMA (pHMM)"] = '|'.join(domain_to_prediction[domain].phmm)
        domain_to_method_to_prediction[domain]["SANDPUMA (prediCAT monophyly)"] = '|'.join(domain_to_prediction[domain].predicat_mp)
        domain_to_method_to_prediction[domain]["SANDPUMA (prediCAT SNN)"] = '|'.join(domain_to_prediction[domain].predicat_snn)


def add_nrpspredictor_data(domain_to_method_to_prediction):
    domain_to_prediction = parse_nrpspredictor_data()
    for domain in domain_to_method_to_prediction:
        domain_to_method_to_prediction[domain]["NRPSPredictor2"] = domain_to_prediction[domain].prediction


def add_adenpredictor_data(domain_to_method_to_prediction):
    domain_to_prediction = parse_adenpredictor_data()
    for domain in domain_to_method_to_prediction:
        domain_to_method_to_prediction[domain]["AdenPredictor"] = domain_to_prediction[domain].prediction

def run():
    args = parse_arguments()

    get_summary(args.i, args.o, args.d, args.m, sandpuma=args.s, adenpredictor=args.a, nrpspredictor=args.n)

if __name__ == "__main__":
    run()