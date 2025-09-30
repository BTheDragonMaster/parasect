import pandas as pd
import scipy.cluster.hierarchy as hc
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.spatial as sp


def write_matrix(matrix, labels, out_file):
    with open(out_file, 'w') as out:
        out.write('Predicted\\Actual\t')
        out.write('\t'.join(labels))
        out.write("\n")
        for i, aa in enumerate(labels):
            out.write(f"{aa}")
            for j, aa_2 in enumerate(labels):
                out.write(f"\t{matrix[i][j]}")
            out.write("\n")


def plot_matrix(matrix, labels, out_file, cluster=True):
    data = pd.DataFrame(matrix, index=labels, columns=labels)

    cmap = sns.color_palette("rocket_r", as_cmap=True)

    if cluster:

        data_correlation = data.T.corr().fillna(0)
        distance_matrix = 1 - data_correlation
        distance_matrix = distance_matrix.clip(lower=0)
        linkage = hc.linkage(sp.distance.squareform(distance_matrix, checks=False), method='average')

        clustermap = sns.clustermap(data, cmap=cmap,
                                    row_linkage=linkage, col_linkage=linkage)
        clustermap.fig.savefig(out_file)
        plt.close(clustermap.fig)
    else:
        ax = sns.heatmap(data, cmap=cmap)
        ax.get_figure().savefig(out_file)
        plt.close(ax.get_figure())
