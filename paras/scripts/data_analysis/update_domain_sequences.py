from sys import argv


from paras.scripts.parsers.fasta import read_fasta
from paras.scripts.parsers.tabular import Tabular


def update_domain_sequences(parasect_dataset, out_file):
    parasect_data = Tabular(parasect_dataset, [0])
    domain_to_seq = {}

    for datapoint in parasect_data.data:
        domain_id = parasect_data.get_value(datapoint, "domain_id")
        domain_to_seq[domain_id] = parasect_data.get_value(datapoint, "sequence")

    domain_ids = list(domain_to_seq.keys())
    domain_ids.sort()
    with open(out_file, 'w') as out:

        for domain_id in domain_ids:
            out.write(f">{domain_id}\n{domain_to_seq[domain_id]}\n")


if __name__ == "__main__":
    update_domain_sequences(argv[1], argv[2])