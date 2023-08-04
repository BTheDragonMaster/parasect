from paras.scripts.parsers.parsers import parse_specificities, parse_substrate_list, parse_domain_list
from pprint import pprint
from sys import argv
import os


def get_cluster_distributions(cluster_file, specificities_file, specificities, out_dir):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    cluster_to_spec_to_count = {}
    domain_to_spec = parse_specificities(specificities_file)
    all_specs = parse_substrate_list(specificities)
    domain_to_cluster = {}

    with open(cluster_file, 'r') as cluster_info:
        for line in cluster_info:
            line = line.strip()
            domain, cluster = line.split('\t')
            cluster = int(cluster)
            if cluster not in cluster_to_spec_to_count:
                cluster_to_spec_to_count[cluster] = {}
                for spec in all_specs:
                    cluster_to_spec_to_count[cluster][spec] = 0

            specs = domain_to_spec[domain]
            count_domain = False
            for spec in specs:
                if spec in all_specs:
                    cluster_to_spec_to_count[cluster][spec] += 1
                    count_domain = True

            if count_domain:
                domain_to_cluster[domain] = cluster

    for spec in all_specs:
        clusters = []
        for cluster, spec_to_count in cluster_to_spec_to_count.items():
            count = spec_to_count[spec]
            if count != 0:
                clusters.append(cluster)
        if len(clusters) < 3:
            print(spec, clusters)

    new_clustering = {}

    for cluster, spec_to_count in cluster_to_spec_to_count.items():
        if cluster in [0, 1]:
            index = 0
        else:
            index = 1

        if index not in new_clustering:
            new_clustering[index] = spec_to_count
        else:
            for spec, count in spec_to_count.items():
                new_clustering[index][spec] += count

    pprint(new_clustering)
    cluster_1 = 0
    cluster_2 = 0
    for spec in all_specs:
        print(f"{spec}\t{new_clustering[0][spec]}\t{new_clustering[1][spec]}")
        cluster_1 += new_clustering[0][spec]
        cluster_2 += new_clustering[1][spec]

    print(f"{cluster_1} datapoints in cluster 1")
    print(f"{cluster_2} datapoints in cluster 2")

    out_cluster = os.path.join(out_dir, "cluster_distribution.txt")

    out_train = os.path.join(out_dir, "train.txt")
    out_test = os.path.join(out_dir, "test.txt")

    write_spec_split(new_clustering, all_specs, out_cluster)
    write_clustering(domain_to_cluster, out_train, out_test)

    return cluster_to_spec_to_count


def write_spec_split(clustering, all_specs, out_file):
    with open(out_file, 'w') as out:
        for spec in all_specs:
            out.write(f"{spec}\t{clustering[0][spec]}\t{clustering[1][spec]}\n")


def write_clustering(domain_to_cluster, out_train, out_test):
    with open(out_train, 'w') as train:
        with open(out_test, 'w') as test:
            for domain, cluster in domain_to_cluster.items():
                if cluster in [0, 1]:
                    train.write(f"{domain}\n")
                else:
                    test.write(f"{domain}\n")




if __name__ == "__main__":
    distribution = get_cluster_distributions(argv[1], argv[2], argv[3], argv[4])
    # pprint(distribution)
