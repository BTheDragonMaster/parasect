from sys import argv

from paras.scripts.parsers.fasta import read_fasta


full_sequences = read_fasta(argv[1])
domain_sequences = read_fasta(argv[2])
mismatches_out = argv[3]
missing_out = argv[4]

uniprot_ids = {}

for seq_id, sequence in full_sequences.items():
    uniprot_id = '.'.join(seq_id.split('.')[:-1])
    uniprot_ids[uniprot_id] = sequence

with open(mismatches_out, 'w') as mismatches:
    with open(missing_out, 'w') as missing:

        for domain, sequence in domain_sequences.items():
            domain_names = domain.split('|')
            domain_found = False
            mismatch_found = False
            mismatch_corrected = False
            mismatching_proteins = []
            mismatching_sequences = []

            for domain_name in domain_names:
                full_sequence = None
                if domain_found:
                    break
                protein_name = '.'.join(domain_name.split('.')[:-1])
                if protein_name in full_sequences:
                    full_sequence = full_sequences[protein_name]

                elif protein_name in uniprot_ids:
                    full_sequence = uniprot_ids[protein_name]

                if full_sequence is not None:
                    if sequence.upper() not in full_sequence.upper():
                        mismatch_found = True
                        mismatching_proteins.append(protein_name)
                        mismatching_sequences.append(sequence)

                    else:
                        if mismatch_found:
                            mismatch_corrected = True
                        domain_found = True

            if not domain_found:
                print(f"Didn't find full-length sequence for {domain}")
                missing.write(f"{domain}\n")

            if mismatch_found and not mismatch_corrected:
                print(f"Sequences don't match for domain {domain} and {'|'.join(mismatching_proteins)}")
                print(mismatching_sequences)
                mismatches.write(f"{domain}\t{'|'.join(mismatching_proteins)}\n")

