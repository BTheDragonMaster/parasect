from sys import argv

from paras.scripts.parsers.fasta import read_fasta
from paras.scripts.parsers.parsers import parse_domain_list


def check_missing_data(domain_list, fasta_file):
    missing_domains = []
    domain_to_seq = read_fasta(fasta_file)
    domains = parse_domain_list(domain_list)
    for domain in domains:
        if domain not in domain_to_seq:
            missing_domains.append(domain)

    return missing_domains


if __name__ == "__main__":
    print(check_missing_data(argv[1], argv[2]))
    