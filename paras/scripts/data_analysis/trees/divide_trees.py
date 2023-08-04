from argparse import ArgumentParser
from ete3 import Tree
import os
from dataclasses import dataclass
from typing import List
from paras.scripts.parsers.parsers import parse_specificities, parse_substrate_list
from paras.scripts.data_processing.divide_into_substrate_groups import SUBSTRATE_GROUPS


def get_group(substrate):
    for group, substrates in SUBSTRATE_GROUPS.items():
        if substrate in substrates:
            return group
    return None


def parse_arguments():
    parser = ArgumentParser()
    parser.add_argument('-n', type=str, required=True, help='Newick tree.')
    parser.add_argument('-o', type=str, required=True, help="Output directory.")
    parser.add_argument('-s', type=str, required=True, help="Specificity file")
    parser.add_argument('-i', type=str, required=True, help="Included substrates")
    parser.add_argument('-g', action='store_true', help="Clade on substrate group")

    args = parser.parse_args()
    return args


def run():

    args = parse_arguments()
    if not os.path.exists(args.o):
        os.mkdir(args.o)
    domain_to_spec = parse_specificities(args.s)
    nwk = args.n
    included_substrates = parse_substrate_list(args.i)
    get_monophyletic_nodes(nwk, domain_to_spec, included_substrates, args.o, args.g)


@dataclass
class Clade:
    domains: List[str]


def get_monophyletic_nodes(nwk, domain_to_spec, included_substrates, out_folder, group):

    with open(nwk, 'r') as tree_file:
        tree_string = tree_file.read()
        tree = Tree(tree_string)

        if not group:

            for l in tree.iter_leaves():
                specificities = domain_to_spec[l.name]
                new_specificities = []
                for specificity in specificities:
                    if specificity in included_substrates:
                        new_specificities.append(specificity)

                if new_specificities:
                    specificity = new_specificities[0]
                else:
                    specificity = specificities[0]

                l.add_feature("specificity", specificity)

            for substrate in included_substrates:
                clades = []
                for node in tree.get_monophyletic(values=[substrate], target_attr="specificity"):
                    clade = Clade([])
                    for descendant in node.iter_leaves():
                        clade.domains.append(descendant.name)
                    clades.append(clade)

                out_file = os.path.join(out_folder, f"{substrate}.txt")

                with open(out_file, 'w') as out:
                    out.write(f"domain\tclade\n")
                    for i, clade in enumerate(clades):
                        for domain in clade.domains:
                            out.write(f"{domain}\t{i}\n")
        else:
            for l in tree.iter_leaves():
                specificities = domain_to_spec[l.name]
                new_specificities = []
                for specificity in specificities:
                    if specificity in included_substrates:
                        new_specificities.append(specificity)

                if new_specificities:
                    specificity = get_group(new_specificities[0])
                else:
                    print(specificity)
                    specificity = specificities[0]

                l.add_feature("specificity", specificity)

            for group, substrates in SUBSTRATE_GROUPS.items():
                print(group)
                clades = []
                for node in tree.get_monophyletic(values=[group], target_attr="specificity"):
                    # print(node.get_ascii(attributes=["specificity", "name"], show_internal=False))

                    clade = Clade([])
                    for descendant in node.iter_leaves():
                        clade.domains.append(descendant.name)
                    clades.append(clade)

                out_file = os.path.join(out_folder, f"{group}.txt")

                written_domains = []

                with open(out_file, 'w') as out:
                    out.write(f"domain\tclade\n")
                    for i, clade in enumerate(clades):
                        for domain in clade.domains:
                            if domain in written_domains:
                                print(domain)
                            out.write(f"{domain}\t{i}\n")
                            written_domains.append(domain)


if __name__ == "__main__":
    run()

