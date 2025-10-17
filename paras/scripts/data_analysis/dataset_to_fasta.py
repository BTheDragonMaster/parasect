from paras.scripts.parsers.tabular import Tabular
from paras.scripts.parsers.fasta import write_fasta
from sys import argv


def data_to_fasta(parasect_dataset, out_file):
    id_to_seq = {}
    parasect_data = Tabular(parasect_dataset, [0])
    for datapoint in parasect_data.data:
        domain_id = parasect_data.get_value(datapoint, "domain_id")
        sequence = parasect_data.get_value(datapoint, "sequence")
        id_to_seq[domain_id] = sequence

    write_fasta(id_to_seq, out_file)


if __name__ == "__main__":
    data_to_fasta(argv[1], argv[2])

