import os
from sys import argv

from paras.scripts.parsers.fasta import read_fasta
from paras.scripts.parsers.parsers import parse_taxonomy
import paras.data.sequence_data.sequences
import paras.data

DOMAINS = os.path.join(os.path.dirname(paras.data.sequence_data.sequences.__file__), "domains.fasta")
TAXONOMY = os.path.join(os.path.dirname(paras.data.__file__), "taxonomy.txt")


def assign_taxonomy(fasta_file, taxonomy_file, out_file):
    id_to_seq = read_fasta(fasta_file)
    protein_to_taxonomy = parse_taxonomy(taxonomy_file)

    no_tax = 0
    dubious_tax = 0
    fungal = 0
    bacterial = 0
    archaea = 0

    domain_to_tax = {}
    print(len(id_to_seq.keys()))

    for seq_id_str in id_to_seq.keys():
        seq_id_list = seq_id_str.split('|')
        protein_list = []
        for seq_id in seq_id_list:
            protein_list.append('.'.join(seq_id.split('.')[:-1]))

        taxonomy_list = []
        for protein in protein_list:
            if protein in protein_to_taxonomy and protein_to_taxonomy[protein] not in taxonomy_list and protein_to_taxonomy[protein]:
                taxonomy_list.append(protein_to_taxonomy[protein])

        if len(taxonomy_list) == 0:
            no_tax += 1
        if len(taxonomy_list) > 1:
            dubious_tax += 1

        if taxonomy_list:
            domain_to_tax[seq_id_str] = taxonomy_list[0]
            if taxonomy_list[0][0] == 'Bacteria':
                bacterial += 1
            elif taxonomy_list[0][0] == "Eukaryota":
                fungal += 1
            else:
                archaea += 1
                print(taxonomy_list)

    with open(out_file, 'w') as out:
        for domain, tax in domain_to_tax.items():
            line = f"{domain}\t" + "\t".join(tax)
            out.write(f"{line}\n")
    print(f"No taxonomy found for {no_tax} domains.")
    print(f"Multiple taxonomies found for {dubious_tax} domains.")
    print(f"{fungal} eukaryotic domains found.")
    print(f"{bacterial} bacterial domains found.")
    print(f"{archaea} unclassified domains found.")


if __name__ == "__main__":
    assign_taxonomy(DOMAINS, TAXONOMY, argv[1])