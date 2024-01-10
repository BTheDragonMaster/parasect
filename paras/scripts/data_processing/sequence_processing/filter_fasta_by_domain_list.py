from sys import argv
import os
from paras.scripts.parsers.fasta import read_fasta
from paras.scripts.parsers.parsers import parse_domain_list
from paras.scripts.parsers.iterate_over_dir import iterate_over_dir


def filter_fasta(alignment_fasta, domain_list, out_path):
    domains_to_keep = parse_domain_list(domain_list)
    domains = []
    sequences = []
    domain_to_seq = read_fasta(alignment_fasta)
    for domain, seq in domain_to_seq.items():
        if domain in domains_to_keep:
            domains.append(domain)
            sequences.append(seq)

    with open(out_path, 'w') as out:
        for i, domain in enumerate(domains):
            out.write(f">{domain}\n{sequences[i]}\n")


def filter_domains(alignment_fasta, domain_list, out_path):
    domains = parse_domain_list(domain_list)
    domain_to_seq = read_fasta(alignment_fasta)
    domains_to_keep = []
    for domain in domains:
        if domain in domain_to_seq:
            domains_to_keep.append(domain)

    with open(out_path, 'w') as out:
        for domain in domains_to_keep:
            out.write(f"{domain}\n")


def filter_domain_folder(alignment_fasta, domain_dir, out_dir):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    for file_name, domain_file in iterate_over_dir(domain_dir, extension='.txt'):
        out_path = os.path.join(out_dir, f"{file_name}.txt")
        filter_domains(alignment_fasta, domain_file, out_path)


if __name__ == "__main__":
    # filter_fasta(argv[1], argv[2], argv[3])
    filter_domains(argv[1], argv[2], argv[3])
    # filter_domain_folder(argv[1], argv[2], argv[3])