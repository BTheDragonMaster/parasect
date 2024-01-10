from ete3 import Tree, AttrFace, TreeStyle
from collections import OrderedDict
from argparse import ArgumentParser
import os

from paras.scripts.parsers.iterate_over_dir import iterate_over_dir
from paras.scripts.parsers.tabular import Tabular
from paras.scripts.parsers.fasta import read_fasta


def parse_arguments():
    parser = ArgumentParser()
    parser.add_argument('-i', type=str, required=True, help='Input directory of nwk files.')
    parser.add_argument('-f', type=str, default=None, help='Input directory of fasta files.')
    parser.add_argument('-s', type=str, required=True, help="Specificity file")
    parser.add_argument('-o', type=str, required=True, help="Output directory.")

    args = parser.parse_args()
    return args


def parse_specificities(data_file):
    domain_to_spec = OrderedDict()
    parasect_data = Tabular(data_file, [0])
    for data_id in parasect_data.data:
        domain_id = parasect_data.get_value(data_id, "domain_id")
        domain_to_spec[domain_id] = parasect_data.get_value(data_id, "specificity").split('|')
    return domain_to_spec


def draw_tree(nwk, nwk_out, specificities, id_to_seq):
    specificities = parse_specificities(specificities)
    with open(nwk, 'r') as tree_file:
        tree_string = tree_file.read()
        tree = Tree(tree_string)

        for l in tree.iter_leaves():
            if id_to_seq is not None:
                l.name = id_to_seq[l.name] + '_' +  '|'.join(specificities[l.name]) + f'_{l.name}'
            else:
                l.name = '|'.join(specificities[l.name]) + f'_{l.name}'

            # create a new label with a color attribute
            N = AttrFace("name")
            N.text = "no"
            # label margins
            N.margin_top = N.margin_bottom = N.margin_left = 4.0
            # labels aligned to the same level
            l.add_face(N, 1, position="branch-right")

        # init tree style
        ts = TreeStyle()
        # remove default labels
        ts.show_leaf_name = False
        # render image on notebook
        tree.render(nwk_out, tree_style=ts)


def run():
    args = parse_arguments()
    if not os.path.exists(args.o):
        os.mkdir(args.o)
    for spec, file_path in iterate_over_dir(args.i, '.nwk'):
        id_to_seq = None
        if args.f is not None:
            fasta_file = os.path.join(args.f, f"{spec}.fasta")
            print(spec)
            id_to_seq = read_fasta(fasta_file)
        svg_name = os.path.join(args.o, f"{spec}.svg")
        draw_tree(file_path, svg_name, args.s, id_to_seq)


if __name__ == "__main__":

    run()