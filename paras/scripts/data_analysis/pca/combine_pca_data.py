from collections import OrderedDict
from sys import argv

from paras.scripts.parsers.parsers import parse_specificities, yield_pocket_features
from paras.scripts.parsers.iterate_over_dir import iterate_over_dir
from paras.scripts.data_analysis.pca.pca import Pca


def run_and_combine_pcas(input_dir, voxel_data, specificities_file, out_path):

    domains = []
    vectors = []

    print("Loading data..")
    # domain_to_pocket = parse_pocket_features(voxel_data)
    domain_to_spec = parse_specificities(specificities_file)
    for domain, pocket, _ in yield_pocket_features(voxel_data):
        domains.append(domain)
        vectors.append(pocket)

    pc_to_vectors = OrderedDict()

    print("Running PCA..")

    categories = []
    for pca_name, model_path in iterate_over_dir(input_dir, '.model'):
        group_name = pca_name.split('pca_')[1]

        pca = Pca.from_model(model_path)

        pca.apply(vectors)
        pc_to_vectors[group_name] = pca.data.vectors
        for i in range(1, pca.pca.n_components + 1):
            categories.append(f"{group_name}_pc{i}")

    with open(out_path, 'w') as out:
        out.write("domain\tspecificity\t")
        out.write('\t'.join(categories))
        out.write('\n')
        for i, domain in enumerate(domains):
            out.write(f"{domain}\t{'|'.join(domain_to_spec[domain])}\t")
            full_vector = []
            for group, vectors in pc_to_vectors.items():
                full_vector += list(vectors[i])
            assert len(categories) == len(full_vector)
            out.write('\t'.join(list(map(str, full_vector))))
            out.write('\n')


if __name__ == "__main__":
    run_and_combine_pcas(argv[1], argv[2], argv[3], argv[4])
