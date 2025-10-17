from paras.scripts.parsers.parsers import parse_domain_list
from paras.scripts.parsers.fasta import read_fasta
from sys import argv


def get_missing(fasta_file, domain_list, out_fasta):
    id_to_seq = read_fasta(fasta_file)
    domains = parse_domain_list(domain_list)
    with open(out_fasta, 'w') as out:
        for domain in domains:
            out.write(f">{domain}\n{id_to_seq[domain]}\n")


if __name__ == "__main__":
    get_missing(argv[1], argv[2], argv[3])
