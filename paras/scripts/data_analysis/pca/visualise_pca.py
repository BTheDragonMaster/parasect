import argparse
import os

from paras.scripts.data_analysis.pca.pca import PcaData


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', type=str, required=True, help="File containing precomputed PCs.")
    parser.add_argument('-v', type=str, default='2d', help="Projection in PCA visualisation")
    parser.add_argument('-t', type=int, default=0, help="Minimal number of occurrences of a substrate for visualising")
    parser.add_argument('-x', type=int, default=1, help="PCA index to visualise in the x-direction")
    parser.add_argument('-y', type=int, default=2, help="PCA index to visualise in the y-direction")
    parser.add_argument('-z', type=int, default=3, help="PCA index to visualise in the z-direction")
    parser.add_argument('-m', type=str, default='min', help="Specifies if number of occurrences is upper or lower bound")
    parser.add_argument('-sv', type=str, nargs='*', default=None, help="Substrates to visualise")
    parser.add_argument('-o', type=str, required=True, help="Path to output directory")
    parser.add_argument('-p', type=str, default=None, help="Prefix")
    parser.add_argument('-d', type=str, nargs='*', default=None, help="Names of domains of interest")
    parser.add_argument('-l', type=str, nargs='*', default=None, help="Labels of domains of interest")

    args = parser.parse_args()
    return args


def run():
    args = parse_arguments()
    if not os.path.exists(args.o):
        os.mkdir(args.o)

    pca_data = PcaData.from_file(args.f)
    print(pca_data.domains)

    if args.v == '2d':
        if args.p:
            png_name = f'{args.p}_pocket_{args.x}_{args.y}.png'
        else:
            png_name = f'pocket_{args.x}_{args.y}.png'
    elif args.v == '3d':
        if args.p:
            png_name = f'{args.p}_pocket_{args.x}_{args.y}_{args.z}.png'
        else:
            png_name = f'pocket_{args.x}_{args.y}_{args.z}.png'

    else:
        raise ValueError("Projection must be 2d or 3d")

    if args.d:
        assert args.l and len(args.l) == len(args.d)

    png_out = os.path.join(args.o, png_name)
    pca_data.visualise(png_out, args.m, args.t, projection=args.v, x=args.x, y=args.y, z=args.z, substrates=args.sv,
                       domains_of_interest=args.d, domain_labels=args.l)


if __name__ == "__main__":
    run()
