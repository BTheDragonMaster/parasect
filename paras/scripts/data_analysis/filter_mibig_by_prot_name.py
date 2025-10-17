from paras.scripts.parsers.fasta import read_fasta, write_fasta
from sys import argv


def filter_mibig_domains(mibig_domain_file, domain_file, out_fasta):
    mibig_to_seq = read_fasta(mibig_domain_file)
    domain_to_seq = read_fasta(domain_file)
    proteins = set()

    for domain in domain_to_seq.keys():
        domain_names = domain.split('|')
        for domain_name in domain_names:
            proteins.add('.'.join(domain_name.split('.')[:-1]))

    print(proteins)

    filtered_to_seq = {}
    for mibig_id in mibig_to_seq.keys():
        protein_id = mibig_id.split('|')[6]
        print(protein_id)
        if protein_id not in proteins:
            filtered_to_seq[mibig_id] = mibig_to_seq[mibig_id]

    write_fasta(filtered_to_seq, out_fasta)


if __name__ == "__main__":
    filter_mibig_domains(argv[1], argv[2], argv[3])


