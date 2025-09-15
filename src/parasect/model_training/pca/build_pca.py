from argparse import ArgumentParser, Namespace

from parasect.core.pca import esm_to_pca


def parse_arguments() -> Namespace:
    """Parse arguments from command line

    :return: Arguments
    :rtype: Namespace
    """
    parser = ArgumentParser(description="Build PCA from ESM embeddings")
    parser.add_argument("-e", '--esm', required=True, type=str,
                        help="Path to folder containing precomputed ESM embeddings")
    parser.add_argument("-o", '--out', required=True, type=str,
                        help="Path to output directory")
    parser.add_argument("-n", "--n_components", type=int, default=500,
                        help="Number of PCA components")
    parser.add_argument("-s", "--store_pca", action="store_true",
                        help="If given, store PCA model")

    arguments = parser.parse_args()

    return arguments


if __name__ == "__main__":
    args = parse_arguments()
    esm_to_pca(esm_dir=args.esm, out_dir=args.out, n_components=args.n_components, store_pca=args.store_pca)
