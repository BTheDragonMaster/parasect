from Bio import SeqIO
from paras.scripts.parsers.fasta import write_fasta


def proteins_from_genbank(gbk_in: str, fasta_out: str) -> None:
    """
    Convert a GBK file to a fasta file containing protein sequences

    Input:
    gbk_in: str, path to .gbk/.gb file
    fasta_out: str, path to output file

    """
    counter = 0
    id_to_seq = {}
    for record in SeqIO.parse(gbk_in, "genbank"):
        for feature in record.features:
            if feature.type == 'CDS':
                if 'translation' in feature.qualifiers:
                    sequence = feature.qualifiers['translation']

                    if 'protein_id' in feature.qualifiers:
                        seq_id = feature.qualifiers['protein_id'][0]
                    elif 'gene_id' in feature.qualifiers:
                        seq_id = feature.qualifiers['gene_id'][0]
                    elif 'locus_tag' in feature.qualifiers:
                        seq_id = feature.qualifiers['locus_tag'][0]
                    else:
                        counter += 1
                        seq_id = f"gene_{counter}"
                    id_to_seq[seq_id] = sequence

    write_fasta(id_to_seq, fasta_out)
