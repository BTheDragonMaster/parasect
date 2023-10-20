import os

from paras.scripts.parsers.fasta import read_fasta

MAPPING_SUFFIX = 'mapping.txt'
FASTA_SUFFIX = 'renamed_fasta.txt'


def rename_sequences(fasta_file, out_dir):
    """
    Rename sequences to integers before running hmmscan

    Save a mapping file which maps these integers to the original sequence IDs
    """
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    mapping_file = os.path.join(out_dir, MAPPING_SUFFIX)
    new_fasta_file = os.path.join(out_dir, FASTA_SUFFIX)

    id_to_seq = read_fasta(fasta_file)
    counter = 0
    with open(new_fasta_file, 'w') as new_fasta:
        with open(mapping_file, 'w') as mapping:
            for seq_id, seq in id_to_seq.items():
                counter += 1
                mapping.write(f"{counter}\t{seq_id}\n")
                new_fasta.write(f">{counter}\n{seq}\n")

    return mapping_file, new_fasta_file


def parse_mapping(mapping_file):
    new_to_original = {}
    with open(mapping_file, 'r') as mapping:
        for line in mapping:
            line = line.strip()
            line_segments = line.split('\t')
            new = line_segments[0]
            original = '\t'.join(line_segments[1:])
            new_to_original[new] = original

    return new_to_original


def reverse_renaming(adenylation_domains, mapping_file):
    new_to_original = parse_mapping(mapping_file)
    for domain in adenylation_domains:
        domain.protein_name = new_to_original[domain.protein_name]
