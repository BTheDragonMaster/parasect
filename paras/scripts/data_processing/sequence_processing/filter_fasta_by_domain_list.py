from sys import argv
from paras.scripts.parsers.fasta import read_fasta
from paras.scripts.parsers.parsers import parse_domain_list


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


if __name__ == "__main__":
    filter_fasta(argv[1], argv[2], argv[3])