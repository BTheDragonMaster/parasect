import os

import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
from joblib import dump, load
from kneed import KneeLocator

from paras.scripts.parsers.tabular import Tabular
from paras.scripts.parsers.parsers import parse_specificities, parse_domain_list


class PcaData:
    def __init__(self, n_components):
        self.domains = []
        self.labels = []
        self.n_components = n_components
        self.vectors = []

    @classmethod
    def from_file(cls, pca_file):
        pca_data = Tabular(pca_file, [0])
        domains = []
        substrates = []
        pca_vectors = []

        for domain in pca_data.data:
            domains.append(pca_data.get_value(domain, 'domain'))
            substrates.append(pca_data.get_value(domain, 'specificity'))
            pca_vectors.append(list(map(float, pca_data.get_row(domain)[2:])))

        n_components = len(pca_vectors[0])
        pca_data = cls(n_components)
        pca_data.set_labels_from_lists(substrates, domains)
        pca_data.vectors = pca_vectors

        return pca_data

    def set_labels_from_file(self, parasect_data, domain_list):
        domain_to_spec = parse_specificities(parasect_data)
        domain_list = parse_domain_list(domain_list)
        for domain in domain_list:
            self.labels.append('|'.join(domain_to_spec[domain]))
            self.domains.append(domain)

    def set_labels_from_lists(self, substrates, domains):
        self.labels = substrates
        self.domains = domains

    def write_transformed_vectors(self, out_file):
        with open(out_file, 'w') as out:
            out.write("domain\tspecificity")
            for i in range(self.n_components):
                out.write(f"\tpc_{i + 1}")
            out.write('\n')

            for i, domain in enumerate(self.domains):
                out.write(f"{domain}\t{self.labels[i]}")
                for value in self.vectors[i]:
                    out.write(f"\t{value}")
                out.write('\n')

    def visualise(self, out_png, threshold_type='min', threshold=0, projection='2d', x=1, y=2, z=3, substrates=None,
                  domains_of_interest=None, domain_labels=None):

        fig = plt.figure(figsize=(8, 6))
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

                substrate_to_count = {}

                for substrate in set(self.labels):
                    if substrates is not None:
                        if substrate in substrates:
                            substrate_to_xy[substrate] = {'x': [],
                                                          'y': [],
                                                          'z': []}
                    else:

                        substrate_to_count[substrate] = self.labels.count(substrate)
                        substrate_to_xy[substrate] = {'x': [],
                                                      'y': [],
                                                      'z': []}

                labels = set()

                for i, substrate in enumerate(self.labels):
                    if threshold_type == 'min':
                        if substrate_to_count[substrate] > threshold:
                            substrate_to_xy[substrate]['x'].append(self.vectors[i][x])
                            substrate_to_xy[substrate]['y'].append(self.vectors[i][y])
                            substrate_to_xy[substrate]['z'].append(self.vectors[i][z])
                            labels.add(substrate)
                    elif threshold_type == 'max':
                        if substrate_to_count[substrate] <= threshold:
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

            plt.gcf().subplots_adjust(right=0.7)

            if projection == '2d':


                for i, substrate in enumerate(filtered_labels):
                    color = colors[i]
                    plt.scatter(substrate_to_xy[substrate]['x'], substrate_to_xy[substrate]['y'],
                                label=substrate, color=color)
                plt.legend(bbox_to_anchor=(1.05, 1), loc=9)
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
                    fig = plt.figure()
                    ax = fig.add_subplot(projection='3d')
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


class Pca:
    def __init__(self, pca):
        self.pca = pca
        self.data = PcaData(pca.n_components)

    def fit(self, vectors):
        self.pca.fit(vectors)

    def save(self, out_pca):

        dump(self.pca, out_pca)

    def apply(self, vectors):
        self.data.vectors = self.pca.transform(vectors)

    def find_optimal_nr_components(self, out_dir):
        x = list(range(1, self.pca.n_components + 1))
        y = self.pca.explained_variance_ratio_
        kneedle = KneeLocator(x, y, S=1.0, curve="convex", direction="decreasing")

        out_knee = os.path.join(out_dir, f"scree_plot.png")
        kneedle.plot_knee()
        plt.savefig(out_knee)
        plt.clf()
        plt.close()

        return kneedle.knee

    def analyse_importance(self, feature_names):
        matrix = np.transpose(self.pca.components_)
        feature_importances = []

        for feature in matrix:
            feature_importances.append(sum(feature))

        feature_to_importance = {}

        for i, importance in enumerate(feature_importances):
            feature_name = feature_names[i]
            feature_to_importance[feature_name] = importance

        return feature_to_importance

    def analyse_importance_per_component(self, feature_names):
        pc_to_feature_to_importance = {}
        for i, component in enumerate(self.pca.components_):
            pc_to_feature_to_importance[i] = {}
            for j, importance in enumerate(component):
                feature_name = feature_names[j]
                pc_to_feature_to_importance[i][feature_name] = importance

        return pc_to_feature_to_importance

    def write_importance_per_component(self, feature_names, out_dir):
        importance_dir = os.path.join(out_dir, 'feature_importances')
        if not os.path.exists(importance_dir):
            os.mkdir(importance_dir)
        pc_to_feature_to_importance = self.analyse_importance_per_component(feature_names)
        for pc, feature_to_importance in pc_to_feature_to_importance.items():
            out_path = os.path.join(importance_dir, f'pc_{pc}_features.txt')
            with open(out_path, 'w') as out:
                for feature, importance in feature_to_importance.items():
                    out.write(f"{feature}\t{importance}\n")

    @classmethod
    def from_model(cls, pca_file):
        pca = load(pca_file)
        pca = cls(pca)
        return pca
