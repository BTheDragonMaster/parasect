import os
from sys import argv

from paras.scripts.parsers.external_parsers import parse_sandpuma_data
from paras.scripts.parsers.parsers import parse_domain_list, parse_test_results, parse_specificities

import paras.data.structure_data
import paras.data.train_test_splits
import paras.data.benchmarking
import paras.data

benchmarking_dir = os.path.dirname(paras.data.benchmarking.__file__)


SUBSTRATE_FILE = os.path.join(os.path.dirname(paras.data.structure_data.__file__), 'included_substrates.txt')
DOMAINS_CLASS_SANDPUMA = os.path.join(benchmarking_dir, 'benchmark_class_filtered.txt')
DOMAINS_PHYLOGENY_SANDPUMA = os.path.join(benchmarking_dir, 'benchmark_phylogeny_filtered.txt')
DOMAINS_CLASS = os.path.join(os.path.dirname(paras.data.train_test_splits.__file__), 'class/test_filtered.txt')
DOMAINS_PHYLOGENY = os.path.join(os.path.dirname(paras.data.train_test_splits.__file__), 'phylogeny/test_filtered.txt')
PARAS_RESULTS_DIR = os.path.join(benchmarking_dir, 'paras_results')
PARASECT_RESULTS_DIR = os.path.join(benchmarking_dir, 'parasect_results')
PARASECT_DATASET = os.path.join(os.path.dirname(paras.data.__file__), "parasect_dataset.txt")


def update_data(results_dir, method, domain_to_method_to_prediction, prefix="phylogeny"):

    structure_results = os.path.join(results_dir,
                                     f"{prefix}_structure/test_performance/test_results.txt")
    sequence_results = os.path.join(results_dir,
                                    f"{prefix}_sequence/test_performance/test_results.txt")
    sequence_structure_results = os.path.join(results_dir,
                                              f"{prefix}_sequence_structure/test_performance/test_results.txt")

    files = [structure_results,
             sequence_results,
             sequence_structure_results]

    methods = [f"{method} (structure)",
               f"{method} (sequence)",
               f"{method} (sequence/structure)"]

    for i, results_file in enumerate(files):
        method = methods[i]
        domain_to_prediction = parse_test_results(results_file)

        for domain in domain_to_method_to_prediction:

            domain_to_method_to_prediction[domain][method] = domain_to_prediction[domain]


def get_summary(out_file, prefix="phylogeny", sandpuma=True):

    if sandpuma:
        domains_phylogeny = DOMAINS_PHYLOGENY_SANDPUMA
        domains_class = DOMAINS_CLASS_SANDPUMA
    else:
        domains_phylogeny = DOMAINS_PHYLOGENY
        domains_class = DOMAINS_CLASS

    domain_to_spec = parse_specificities(PARASECT_DATASET)

    if prefix == 'phylogeny':
        domain_list = parse_domain_list(domains_phylogeny)
    elif prefix == 'class':
        domain_list = parse_domain_list(domains_class)
    else:
        raise ValueError("Prefix must be either phylogeny or class.")

    domain_to_method_to_prediction = {}
    for domain in domain_list:
        domain_to_method_to_prediction[domain] = {}
        domain_to_method_to_prediction[domain]["substrate"] = '|'.join(domain_to_spec[domain])

    update_data(PARASECT_RESULTS_DIR, "PARASECT", domain_to_method_to_prediction, prefix)
    update_data(PARAS_RESULTS_DIR, "PARAS", domain_to_method_to_prediction, prefix)
    if sandpuma:
        add_sandpuma_data(domain_to_method_to_prediction)

    with open(out_file, 'w') as out:
        methods = []
        for method in domain_to_method_to_prediction[list(domain_to_method_to_prediction.keys())[0]]:
            methods.append(method)

        methods.sort()
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


if __name__ == "__main__":
    out_dir = argv[1]
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    get_summary(os.path.join(out_dir, "phylogeny.txt"), prefix="phylogeny", sandpuma=False)
    get_summary(os.path.join(out_dir, "class.txt"), prefix="class", sandpuma=False)