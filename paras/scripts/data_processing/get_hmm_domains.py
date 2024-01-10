from paras.scripts.parsers.fasta import read_fasta
from paras.scripts.parsers.parsers import parse_domain_list
from sys import argv

def get_hmm_domain_list(full_domain_list, hmm_fasta, out_file):
    domains = parse_domain_list(full_domain_list)
    hmm_domains = set(read_fasta(hmm_fasta).keys())
    with open(out_file, 'w') as out:
        for domain in domains:
            if domain in hmm_domains:
                out.write(f"{domain}\n")


if __name__ == "__main__":
    get_hmm_domain_list(argv[1], argv[2], argv[3])