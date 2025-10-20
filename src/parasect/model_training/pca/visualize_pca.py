import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
import argparse
import os

from sqlalchemy.orm import Session
from sqlalchemy import create_engine

from parasect.core.tabular import Tabular
from parasect.core.parsing import parse_list
from parasect.database.build_database import Substrate
from parasect.database.query_database import get_domains_from_synonym, get_substrates_from_name, \
    get_domains_from_taxonomic_rank
from parasect.model_training.train_test_splits.substrate_selection import SubstrateSelectionMode, \
    map_domains_to_substrates
from parasect.core.taxonomy import Rank


class PcaData:
    def __init__(self, n_components):
        self.domains = []
        self.labels = []
        self.n_components = n_components
        self.vectors = []

    @classmethod
    def from_file(cls, session, pca_file, included_substrates: set[Substrate],
                  selection_mode: SubstrateSelectionMode = SubstrateSelectionMode.FIRST_ONLY,
                  by_clade: bool = False):
        assert selection_mode in [SubstrateSelectionMode.FIRST_ONLY, SubstrateSelectionMode.FIRST_VALID]
        pca_data = Tabular(pca_file)
        domains = []
        substrates = []
        pca_vectors = []

        for domain_name in pca_data.rows:
            domain_synonym = domain_name.split('|')[0]
            domain = get_domains_from_synonym(session, domain_synonym)[0]
            domains.append(domain)
            pca_vectors.append(list(map(float, pca_data.get_row_values(domain_name)[1:])))

        domain_to_substrates = map_domains_to_substrates(domains, included_substrates, selection_mode)

        filtered_domains = []
        filtered_vectors = []

        fungal_domains = get_domains_from_taxonomic_rank(session, Rank.KINGDOM, "Fungi")
        bacterial_domains = get_domains_from_taxonomic_rank(session, Rank.DOMAIN, "Bacteria")

        for i, domain in enumerate(domains):
            if domain_to_substrates[domain] is not None:
                filtered_domains.append(domain)
                if by_clade:
                    if domain in fungal_domains:
                        substrates.append("Fungal")
                    elif domain in bacterial_domains:
                        substrates.append("Bacterial")
                    else:
                        substrates.append("Unknown")
                else:
                    substrates.append(domain_to_substrates[domain].name)
                filtered_vectors.append(pca_vectors[i])

        n_components = len(pca_vectors[0])
        pca_data = cls(n_components)
        pca_data.set_labels_from_lists(substrates, filtered_domains)
        pca_data.vectors = filtered_vectors

        return pca_data

    def set_labels_from_lists(self, substrates, domains):
        self.labels = substrates
        self.domains = domains

    def visualise(self, out_png, projection='2d', x=1, y=2, z=3, substrates=None,
                  domains_of_interest=None, domain_labels=None):

        fig = plt.figure(figsize=(15, 8))
        x -= 1
        y -= 1
        z -= 1

        # plt.rcParams["figure.figsize"] = (10, 8)

        if not domains_of_interest:

            substrate_to_xy = {}

            if substrates:
                for substrate in set(self.labels):
                    if substrate in substrates:
                        substrate_to_xy[substrate] = {'x': [],
                                                      'y': [],
                                                      'z': []}
                labels = set()
                for i, substrate in enumerate(self.labels):
                    if substrate in substrates:
                        substrate_to_xy[substrate]['x'].append(self.vectors[i][x])
                        substrate_to_xy[substrate]['y'].append(self.vectors[i][y])
                        substrate_to_xy[substrate]['z'].append(self.vectors[i][z])
                        labels.add(substrate)

            else:

                for substrate in set(self.labels):
                    if substrates is not None:
                        if substrate in substrates:
                            substrate_to_xy[substrate] = {'x': [],
                                                          'y': [],
                                                          'z': []}
                    else:

                        substrate_to_xy[substrate] = {'x': [],
                                                      'y': [],
                                                      'z': []}

                labels = set()

                for i, substrate in enumerate(self.labels):
                    substrate_to_xy[substrate]['x'].append(self.vectors[i][x])
                    substrate_to_xy[substrate]['y'].append(self.vectors[i][y])
                    substrate_to_xy[substrate]['z'].append(self.vectors[i][z])
                    labels.add(substrate)

            if len(labels) <= 22:
                colors = list(cm.tab20(np.linspace(0, 1, 20)))
                colors = ['black'] + colors + ['blue']
            else:
                colors = cm.nipy_spectral(np.linspace(0, 1, len(labels)))

            filtered_labels = sorted(labels)

            plt.gcf().subplots_adjust(left=0.05, right=0.45)

            if projection == '2d':

                for i, substrate in enumerate(filtered_labels):
                    color = colors[i]
                    plt.scatter(substrate_to_xy[substrate]['x'], substrate_to_xy[substrate]['y'],
                                label=substrate, color=color)
                plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
                plt.savefig(out_png)

            elif projection == '3d':
                ax = fig.add_subplot(projection='3d')
                for i, substrate in enumerate(filtered_labels):
                    color = colors[i]
                    ax.scatter(substrate_to_xy[substrate]['x'],
                               substrate_to_xy[substrate]['y'],
                               substrate_to_xy[substrate]['z'],
                               label=substrate, color=color)

                ax.set_xlabel(f'PC {x + 1}')
                ax.set_ylabel(f'PC {y + 1}')
                ax.set_zlabel(f'PC {z + 1}')

                ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
                fig.savefig(out_png)

            else:
                raise ValueError("Projection must be 2d or 3d")

        else:
            if domain_labels:
                assert len(domains_of_interest) == len(domain_labels)

            if len(domains_of_interest) == 2:
                colors = ["lightsteelblue", "darksalmon"] #["gainsboro"]
            elif len(domains_of_interest) <= 20:
                colors = list(cm.tab20(np.linspace(0, 1, 20)))
            else:
                colors = cm.nipy_spectral(np.linspace(0, 1, len(domains_of_interest)))

            xs = []
            ys = []
            zs = []

            for vector in self.vectors:
                xs.append(vector[x])
                ys.append(vector[y])
                zs.append(vector[z])
            if projection == '2d':
                plt.scatter(xs, ys, color="gainsboro")
            elif projection == '3d':
                fig = plt.figure()
                ax = fig.add_subplot(projection='3d')
                ax.scatter(xs, ys, zs, color="gainsboro")
            else:
                raise ValueError("Projection must be 2d or 3d")

            for i, domain in enumerate(domains_of_interest):
                color = colors[i]
                domain_index = self.domains.index(domain)
                if domain_labels:
                    label = domain_labels[i]
                else:
                    label = domain
                domain_vector = self.vectors[domain_index]
                xs = [domain_vector[x]]
                ys = [domain_vector[y]]
                zs = [domain_vector[z]]

                if projection == '2d':
                    plt.scatter(xs, ys, color=color, label=label)
                elif projection == '3d':
                    ax.scatter(xs, ys, zs, color=color, label=label)
                else:
                    raise ValueError("Projection must be 2d or 3d")

            if projection == '2d':
                legend = plt.legend(loc='upper right')
                plt.savefig(out_png, bbox_extra_artists=(legend, ), bbox_inches='tight')
                plt.xlabel(f'PC {x + 1}')
                plt.ylabel(f'PC {y + 1}')
            elif projection == '3d':
                ax.set_xlabel(f'PC {x + 1}')
                ax.set_ylabel(f'PC {y + 1}')
                ax.set_zlabel(f'PC {z + 1}')

                ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
                fig.savefig(out_png)
            else:
                raise ValueError("Projection must be 2d or 3d")

        plt.show()
        plt.clf()




def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-pca', '--pca_file', type=str, required=True, help="File containing precomputed PCs.")
    parser.add_argument('-v', '--projection', type=str, default='2d', help="Projection in PCA visualisation")
    parser.add_argument('-x', type=int, default=1, help="PCA index to visualise in the x-direction")
    parser.add_argument('-y', type=int, default=2, help="PCA index to visualise in the y-direction")
    parser.add_argument('-z', type=int, default=3, help="PCA index to visualise in the z-direction")
    parser.add_argument('-s', '--shown_substrates', type=str, nargs='*', default=None, help="Substrates to show")
    parser.add_argument('-o', '--out', type=str, required=True, help="Path to output directory")
    parser.add_argument('-p', '--prefix', type=str, default=None, help="Prefix")
    parser.add_argument('-d', '--highlighted_domains', type=str, nargs='*', default=None,
                        help="Domains to highlight")
    parser.add_argument('-l', '--domain_labels', type=str, nargs='*', default=None, help="Labels of domains of interest")
    parser.add_argument('-i', '--included_substrates', type=str, required=True, help="Path to file with included substrates")
    parser.add_argument('-db', '--database', type=str, required=True, help="Path to sql database")
    parser.add_argument('-c', '--visualise_clades', action='store_true',
                        help="If given, visualise taxonomy instead of substrates")

    args = parser.parse_args()
    return args


def run():
    args = parse_arguments()
    if not os.path.exists(args.out):
        os.mkdir(args.out)

    engine = create_engine(f"sqlite:///{args.database}")

    with Session(engine) as session:
        included_substrate_names = parse_list(args.included_substrates)
        included_substrates = set()
        for name in included_substrate_names:
            included_substrates.add(get_substrates_from_name(session, name)[0])

        pca_data = PcaData.from_file(session, args.pca_file, included_substrates, by_clade=args.visualise_clades)

        if args.projection == '2d':
            if args.prefix:
                png_name = f'{args.prefix}_pocket_{args.x}_{args.y}.png'
            else:
                png_name = f'pocket_{args.x}_{args.y}.png'
        elif args.projection == '3d':
            if args.prefix:
                png_name = f'{args.prefix}_pocket_{args.x}_{args.y}_{args.z}.png'
            else:
                png_name = f'pocket_{args.x}_{args.y}_{args.z}.png'

        else:
            raise ValueError("Projection must be 2d or 3d")

        if args.highlighted_domains:
            assert args.domain_labels and len(args.domain_labels) == len(args.highlighted_domains)

        png_out = os.path.join(args.out, png_name)
        pca_data.visualise(png_out, projection=args.projection, x=args.x, y=args.y, z=args.z,
                           substrates=args.shown_substrates,
                           domains_of_interest=args.highlighted_domains, domain_labels=args.domain_labels)


if __name__ == "__main__":
    run()
