from collections import OrderedDict
from sys import argv

from paras.scripts.parsers.parsers import parse_substrate_list
from paras.scripts.parsers.fasta import read_fasta, write_fasta
from paras.scripts.feature_extraction.sequence_feature_extraction.sequence_labels import SEPARATOR_1


def relabel_from_inclusion_limit(domain_list, domain_to_substrates, substrate_counts, limit):
    domain_to_filtered = {}
    for domain in domain_list:
        for substrate in domain_to_substrates[domain]:
            if substrate_counts[substrate] >= limit:
                if domain not in domain_to_filtered:
                    domain_to_filtered[domain] = []
                domain_to_filtered[domain].append(substrate)

    return domain_to_filtered


def relabel_from_substrate_list(domain_list, domain_to_substrates, substrate_file):
    included_substrates = parse_substrate_list(substrate_file)

    domain_to_filtered = {}
    for domain in domain_list:
        for substrate in domain_to_substrates[domain]:
            if substrate in included_substrates:
                if domain not in domain_to_filtered:
                    domain_to_filtered[domain] = []
                domain_to_filtered[domain].append(substrate)

    return domain_to_filtered


def strip_domain_info(paras_fasta, out_file, separator=SEPARATOR_1):
    new_id_to_seq = OrderedDict()
    old_id_to_seq = read_fasta(paras_fasta)
    for old_id, seq in old_id_to_seq.items():
        new_id = '|'.join(old_id.split(separator)[:-2])
        new_id_to_seq[new_id] = seq

    write_fasta(new_id_to_seq, out_file)


if __name__ == "__main__":
    strip_domain_info(argv[1], argv[2])
