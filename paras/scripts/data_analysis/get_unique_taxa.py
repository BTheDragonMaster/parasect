from paras.scripts.parsers.parsers import parse_taxonomy, parse_domain_list
from sys import argv
from pprint import pprint
import os

def get_taxa_counts(domain_file, taxonomy_file):
    domains = parse_domain_list(domain_file)
    protein_to_tax = parse_taxonomy(taxonomy_file)

    tax_to_count = {}
    for domain in domains:
        protein_string = '.'.join(domain.split('.')[:-1])
        proteins = protein_string.split('|')
        for protein in proteins:
            if protein in protein_to_tax:

                tax_string = '|'.join(protein_to_tax[protein][:3])
                tax_string = tax_string.replace('Bacillota', 'Firmicutes')
                tax_string = tax_string.replace('Proteobacteria', 'Pseudomonadota')
                tax_string = tax_string.replace('Actinobacteria', 'Actinomycetota|Actinomycetes')
                tax_string = tax_string.replace('Cyanobacteria', 'Cyanobacteriota')
                tax_string = '|'.join(tax_string.split('|')[:3])
                if tax_string not in tax_to_count:
                    tax_to_count[tax_string] = 0
                tax_to_count[tax_string] += 1
                break

    return tax_to_count


class Node:
    def __init__(self, parent=None):
        self.parent = parent
        self.children = []
        self.nr_domains = 0

def iterate_over_tree(tree, start_node, parent=None):
    yield parent, start_node

    for child in tree[start_node]:
        yield from iterate_over_tree(tree[start_node], child, start_node)


def build_tree(domain_file, taxonomy_file, out_folder):
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)
    tax_to_count = get_taxa_counts(domain_file, taxonomy_file)
    max_length = 0
    for tax in tax_to_count.keys():
        tax_length = len(tax.split('|'))
        if tax_length > max_length:
            max_length = tax_length

    tree = {'all': {}}
    level_to_count = {'all': 0}
    node_to_root = {}

    for i in range(max_length):
        for tax_string, count in tax_to_count.items():
            tax = tax_string.split('|')
            if i == 0:
                level_to_count['all'] += count

            try:
                level = tax[i]
                node_to_root[level] = tax[0]

                if level not in level_to_count:
                    level_to_count[level] = 0
                level_to_count[level] += count
                current_tax = tree['all']
                for j in range(i):
                    current_tax = current_tax[tax[j]]

                if level not in current_tax:
                    current_tax[level] = {}


            except IndexError:
                continue

    pprint(level_to_count)
    pprint(tree)

    tree_out = os.path.join(out_folder, 'taxonomy_tree.txt')
    domain_counts = os.path.join(out_folder, 'domain_counts.txt')
    colours_out = os.path.join(out_folder, 'colours.txt')

    with open(tree_out, 'w') as out:
        for parent, node in iterate_over_tree(tree, 'all'):
            out.write(f"{parent}\t{node}\n")
    with open(domain_counts, 'w') as out:
        for level, count in level_to_count.items():
            out.write(f"{level}\t{count}\n")

    with open(colours_out, 'w') as out:
        for node, root in node_to_root.items():
            if root == 'Bacteria':
                colour = 1
            elif root == "Eukaryota":
                colour = 2
            else:
                colour = 3
            out.write(f"{node}\t{colour}\n")


if __name__ == "__main__":
    # get_taxa_counts(argv[1], argv[2])
    build_tree(argv[1], argv[2], argv[3])
