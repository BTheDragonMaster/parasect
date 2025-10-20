import os
from typing import Optional
from sklearn.decomposition import PCA
import numpy as np
from joblib import dump

from parasect.core.parsing import iterate_over_dir, parse_esm_embedding


def esm_to_pca(esm_dir: str, out_dir: str, n_components: int = 500, store_pca: bool = False,
               prefix: str = "all_domains", domain_names: Optional[list[str]] = None) -> None:
    """Store transformed features to output

    :param esm_dir: directory containing files, with each file containing precomputed esm embedding for a domain extended signature
    :type esm_dir: str
    :param out_dir: path to output directory
    :type out_dir: str
    :param n_components: number of principal components to keep. Default: 500
    :type n_components: int
    :param store_pca: if True, store PCA model as well as transformed features. Default: False
    :type store_pca: bool
    :param prefix: prefix for file
    :type prefix: str
    :param domain_names: if given, only include these domains
    :type domain_names: Optional[list[str]]
    """
    domain_nr = 0
    for domain_name, domain_path in iterate_over_dir(esm_dir, "_esm.tsv"):
        if domain_names is None or domain_name in domain_names:
            domain_nr += 1
    feature_nr = 1280 * 34
    esm_array = np.zeros(shape=(domain_nr, feature_nr))

    counter = 0
    domains: list[str] = []
    print("Loading data...")
    for domain_name, domain_path in iterate_over_dir(esm_dir, "_esm.tsv"):
        if domain_names is None or domain_name in domain_names:
            domains.append(domain_name)
            esm_embedding = parse_esm_embedding(domain_path)
            esm_array[counter] = esm_embedding
            counter += 1

    pca = PCA(n_components=n_components, random_state=10012025)
    print("Transforming data...")
    pca_features = pca.fit_transform(esm_array)

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    pcs = os.path.join(out_dir, f"{prefix}_pcs_{n_components}.txt")

    with open(pcs, 'w') as out:
        out.write("domain_name")
        for i in range(n_components):
            out.write(f'\tpc_{i + 1}')
        out.write('\n')

        for i, embedding in enumerate(pca_features):
            domain_name = domains[i]
            out.write(domain_name)
            for feature in embedding:
                out.write(f"\t{feature}")
            out.write('\n')

    if store_pca:
        pca_path = os.path.join(out_dir, f"{n_components}.pca")
        dump(pca, pca_path)
