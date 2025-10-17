from sys import argv

from paras.scripts.parsers.fasta import read_fasta, write_fasta
from paras.scripts.parsers.tabular import Tabular


def remove_missing(parasect_dataset, domain_fasta, new_domain_fasta):
    parasect_data = Tabular(parasect_dataset, [0])
    domain_ids = parasect_data.get_column('domain_id')
    id_to_seq = read_fasta(domain_fasta)
    new_id_to_seq = {}

    for seq_id, seq in id_to_seq.items():
        if seq_id in domain_ids:
            new_id_to_seq[seq_id] = seq
        else:
            print(seq_id)

    print("======")

    for seq_id in domain_ids:
        if seq_id not in new_id_to_seq:
            print(seq_id)
    write_fasta(new_id_to_seq, new_domain_fasta)


if __name__ == "__main__":
    remove_missing(argv[1], argv[2], argv[3])