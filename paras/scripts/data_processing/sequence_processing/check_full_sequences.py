from sys import argv

from paras.scripts.parsers.fasta import read_fasta


full_sequences = read_fasta(argv[1])
domain_sequences = read_fasta(argv[2])

for domain, sequence in domain_sequences.items():
    domain_names = domain.split('|')
    for domain_name in domain_names:
        protein_name = '.'.join(domain_name.split('.')[:-1])
        if protein_name in full_sequences:
            full_sequence = full_sequences[protein_name]
            if sequence.upper() not in full_sequence.upper():
                print(f"Sequences don't match for domain {domain} and {protein_name}")
                print(sequence)

        else:
            print(f"Didn't find full-length sequence for {protein_name}")

