from sys import argv

from paras.scripts.parsers.fasta import read_fasta
from paras.scripts.parsers.tabular import Tabular


def inventorise_sequences(sequence_files, specificity_files, file_labels, out_file):
    id_to_domain = {}
    for i, sequence_file in enumerate(sequence_files):
        label = file_labels[i]
        print(label)
        id_to_seq = read_fasta(sequence_file)
        spec_file = specificity_files[i]
        specs = Tabular(spec_file, [0])
        print(specs.categories)
        for seq_id, sequence in id_to_seq.items():
            specificities = specs.get_value(seq_id, "specificity")
            if seq_id not in id_to_domain:
                id_to_domain[seq_id] = Domain(seq_id)
            id_to_domain[seq_id].add_sequence(label, sequence)
            id_to_domain[seq_id].add_specificity(label, specificities)

    with open(out_file, 'w') as out:
        out.write("domain_id\tmibig_seq\tmibig_specificity\tparas_seq\tparas_specificity\tsandpuma_seq\tsandpuma_specificity\tnrpspredictor_seq\tnrpspredictor_specificity\n")
        for domain_id, domain in id_to_domain.items():
            out.write(f"{domain_id}\t{domain.mibig_seq}\t{domain.mibig_spec}\t{domain.paras_seq}\t{domain.paras_spec}\t{domain.sandpuma_seq}\t{domain.sandpuma_spec}\t{domain.nrpspredictor_seq}\t{domain.nrpspredictor_spec}\n")


class Domain:

    def __init__(self, seq_id):
        self.seq_id = seq_id
        self.mibig_seq = ''
        self.paras_seq = ''
        self.sandpuma_seq = ''
        self.nrpspredictor_seq = ''

        self.mibig_spec = ''
        self.paras_spec = ''
        self.sandpuma_spec = ''
        self.nrpspredictor_spec = ''

    def add_sequence(self, label, sequence):
        if label == 'mibig':
            self.mibig_seq = sequence
        elif label == 'paras':
            self.paras_seq = sequence
        elif label == 'sandpuma':
            self.sandpuma_seq = sequence
        elif label == 'nrpspredictor':
            self.nrpspredictor_seq = sequence

    def add_specificity(self, label, specificity):
        if label == 'mibig':
            self.mibig_spec = specificity
        elif label == 'paras':
            self.paras_spec = specificity
        elif label == 'sandpuma':
            self.sandpuma_spec = specificity
        elif label == 'nrpspredictor':
            self.nrpspredictor_spec = specificity


if __name__ == "__main__":
    sequence_files = argv[1:5]
    specificity_files = argv[5:9]
    out_file = argv[9]
    file_labels = ["mibig", "paras", "sandpuma", "nrpspredictor"]
    inventorise_sequences(sequence_files, specificity_files, file_labels, out_file)
