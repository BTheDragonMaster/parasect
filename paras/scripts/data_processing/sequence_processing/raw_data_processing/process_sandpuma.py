import os
from sys import argv

from paras.scripts.parsers.fasta import read_fasta


def process_sandpuma(fasta, out_dir):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    fasta_out_file = os.path.join(out_dir, 'domains.fasta')
    spec_out_file = os.path.join(out_dir, 'specificities.txt')
    id_to_seq = read_fasta(fasta)
    with open(spec_out_file, 'w') as spec_out:
        with open(fasta_out_file, 'w') as fasta_out:
            for seq_id, seq in id_to_seq.items():
                domain_id, specificity, protein_name, adomain_nr = seq_id.split('\t')[:4]
                identifier = f"{protein_name}.{adomain_nr}"
                fasta_out.write(f">{identifier}\n{seq}\n")
                spec_out.write(f"{identifier}\t{specificity}\n")


if __name__ == "__main__":
    process_sandpuma(argv[1], argv[2])