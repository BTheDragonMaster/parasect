from argparse import ArgumentParser, Namespace
import os

from sqlalchemy import create_engine
from sqlalchemy.orm import Session

from parasect.core.pca import esm_to_pca
from parasect.model_training.train_test_splits.domain_scope import DomainScope


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
    parser.add_argument("-db", "--database", type=str, default=None, help="Path to PARASECT database")
    parser.add_argument("-f", "--fungal", action="store_true",
                        help="If given, also create a pca for fungal domains as well. -db must be given")
    parser.add_argument("-b", "--bacterial", action="store_true",
                        help="If given, also create a pca for bacterial domains as well. -db must be given")

    arguments = parser.parse_args()

    if arguments.fungal or arguments.bacterial:
        assert arguments.database is not None

    return arguments


if __name__ == "__main__":
    args = parse_arguments()
    if not os.path.exists(args.out):
        os.mkdir(args.out)
    esm_to_pca(esm_dir=args.esm, out_dir=args.out, n_components=args.n_components, store_pca=args.store_pca,
               prefix="all_domains")

    if args.database:
        engine = create_engine(f"sqlite:///{args.database}")

        with Session(engine) as session:

            if args.fungal:
                fungal_domains_names = [f.get_name() for f in DomainScope.get_domains(session, DomainScope.FUNGAL_ONLY)]
                esm_to_pca(esm_dir=args.esm, out_dir=args.out, n_components=args.n_components, store_pca=args.store_pca,
                           prefix="fungal_domains", domain_names=fungal_domains_names)

            if args.bacterial:
                bacterial_domains = [b.get_name() for b in DomainScope.get_domains(session, DomainScope.BACTERIAL_ONLY)]
                esm_to_pca(esm_dir=args.esm, out_dir=args.out, n_components=args.n_components, store_pca=args.store_pca,
                           prefix="bacterial_domains", domain_names=bacterial_domains)
