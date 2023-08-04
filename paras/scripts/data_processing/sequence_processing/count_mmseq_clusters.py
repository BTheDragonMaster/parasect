from paras.scripts.parsers.parsers import parse_mmseq_clusters, parse_specificities, parse_substrate_list
import argparse
from pprint import pprint


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', type=str, required=True, help="Path to mmseq tsv")
    parser.add_argument('-s', type=str, required=True, help="Path to specificities file")
    parser.add_argument('-i', type=str, required=True, help="Path to included substrates")
    args = parser.parse_args()

    return args


def run():
    args = parse_arguments()

    domain_to_spec = parse_specificities(args.s)
    included_substrates = parse_substrate_list(args.i)
    get_mmseq_metrics(args.m, domain_to_spec, included_substrates)


def get_mmseq_metrics(mmseqs_file, domain_to_spec, included_substrates):
    cluster_to_domain = parse_mmseq_clusters(mmseqs_file)
    cluster_to_spec_counts = {}
    for cluster in cluster_to_domain:
        cluster_to_spec_counts[cluster] = {}
        for substrate in included_substrates:
            cluster_to_spec_counts[cluster][substrate] = 0

    for cluster, domains in cluster_to_domain.items():
        for domain in domains:
            substrates = domain_to_spec[domain]
            for substrate in substrates:
                if substrate in included_substrates:
                    cluster_to_spec_counts[cluster][substrate] += 1

    pprint(cluster_to_spec_counts)

    print("Number of clusters:", len(cluster_to_domain))


if __name__ == "__main__":
    run()
