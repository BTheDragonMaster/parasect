from paras.scripts.parsers.tabular import Tabular
from sys import argv

def table_to_fasta(parasect_data, fasta_out):
    with open(fasta_out, 'w') as out:
        parasect_data = Tabular(parasect_data, [0])
        for prot_id in parasect_data.data:
            out.write(f">{parasect_data.get_value(prot_id, 'domain_id')}\n{parasect_data.get_value(prot_id, 'sequence')}\n")

def table_to_fasta_substrate(parasect_data, fasta_out, substrate):
    with open(fasta_out, 'w') as out:
        parasect_data = Tabular(parasect_data, [0])
        for prot_id in parasect_data.data:
            substrates = parasect_data.get_value(prot_id, 'specificity').split('|')
            if substrate in substrates:
                out.write(f">{parasect_data.get_value(prot_id, 'domain_id')}\n{parasect_data.get_value(prot_id, 'sequence')}\n")


if __name__ == "__main__":
    #table_to_fasta(argv[1], argv[2])
    table_to_fasta_substrate(argv[1], argv[2], argv[3])

