from paras.scripts.parsers.parsers import parse_specificities
from paras.scripts.parsers.fasta import read_fasta

from sys import argv


def rename_nrpspredictor_domains(fasta_file, parasect_dataset, out_file):
    domain_to_specificities = parse_specificities(parasect_dataset)
    domain_to_full_name = {}

    for domain in domain_to_specificities:
        for domain_name in domain.split('|'):
            domain_to_full_name[domain_name] = domain
    id_to_seq = read_fasta(fasta_file)
    with open(out_file, 'w') as out:
        for seq_id, seq in id_to_seq.items():
            if seq_id in domain_to_full_name:
                full_name = domain_to_full_name[seq_id]
                out.write(f">{full_name}\n{seq}\n")
            else:
                print(f"Couldn't find full name for domain {seq_id}")


def make_nrpspredictor_domain_list(fasta_file, parasect_dataset, out_file):
    domain_to_specificities = parse_specificities(parasect_dataset)
    domain_to_full_name = {}

    for domain in domain_to_specificities:
        for domain_name in domain.split('|'):
            domain_to_full_name[domain_name] = domain
    id_to_seq = read_fasta(fasta_file)
    with open(out_file, 'w') as out:
        for seq_id, seq in id_to_seq.items():
            if seq_id in domain_to_full_name:
                full_name = domain_to_full_name[seq_id]
                out.write(f"{full_name}\n")
            else:
                print(f"Couldn't find full name for domain {seq_id}")


if __name__ == "__main__":
    make_nrpspredictor_domain_list(argv[1], argv[2], argv[3])

