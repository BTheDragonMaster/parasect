import os
from sys import argv

from paras.scripts.parsers.fasta import read_fasta, write_fasta
from paras.scripts.parsers.parsers import parse_domain_list
from paras.scripts.parsers.tabular import Tabular


def get_sandpuma_test_set(sandpuma_fasta, parasect_dataset, domain_file, out_dir):
    sandpuma_to_seq = read_fasta(sandpuma_fasta)
    sandpuma_domains = list(sandpuma_to_seq.keys())
    parasect_data = Tabular(parasect_dataset, [0])

    domain_list = parse_domain_list(domain_file)

    sandpuma_test = []
    test_to_sequence = {}

    for domain in domain_list:
        in_sandpuma = False
        for domain_name in domain.split('|'):
            if domain_name in sandpuma_domains:
                in_sandpuma = True

        if not in_sandpuma:
            sandpuma_test.append(domain)
            sequence = parasect_data.get_value((domain,), 'sequence')
            test_to_sequence[domain] = sequence

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    fasta_path = os.path.join(out_dir, "sandpuma_test.fasta")
    domain_path = os.path.join(out_dir, "sandpuma_test.txt")

    write_fasta(test_to_sequence, fasta_path)

    with open(domain_path, 'w') as out:
        for domain in sandpuma_test:
            out.write(f"{domain}\n")


if __name__ == "__main__":
    get_sandpuma_test_set(argv[1], argv[2], argv[3], argv[4])