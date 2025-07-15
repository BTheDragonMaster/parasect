from paras.scripts.parsers.fasta import read_fasta
from Bio import SeqIO, Entrez
from paras.scripts.parsers.parsers import parse_taxonomy
import os
from urllib.error import HTTPError
import time
from sys import argv

def sequences_are_equivalent(seq1: str, seq2: str, min_overlap_ratio: float = 0.8) -> bool:
    """
    Check if two DNA sequences represent the same sequence with partial overlaps.

    Parameters:
    - seq1, seq2: DNA sequences (strings)
    - min_overlap_ratio: minimum fraction of the shorter sequence length to consider as valid overlap

    Returns:
    - True if sequences overlap sufficiently and match in the overlapping region.
    - False otherwise.
    """

    def overlap_length(a, b):
        max_overlap = min(len(a), len(b))
        for length in range(max_overlap, 0, -1):
            if a[-length:] == b[:length]:
                return length
        return 0

    overlap_1_2 = overlap_length(seq1, seq2)
    overlap_2_1 = overlap_length(seq2, seq1)

    min_len = min(len(seq1), len(seq2))
    threshold = int(min_len * min_overlap_ratio)

    if overlap_1_2 >= threshold:
        return True
    if overlap_2_1 >= threshold:
        return True
    if seq1 in seq2:
        return True
    if seq2 in seq1:
        return True

    return False


def chunked_iterable(iterable, size):
    for i in range(0, len(iterable), size):
        yield iterable[i:i + size]


def get_taxid(fasta_file, out_file, chunk_size=200, delay=0.4):
    id_to_seq = read_fasta(fasta_file)
    seq_ids = set()

    for seq_id_str in id_to_seq.keys():
        parts = seq_id_str.split('|')
        if len(parts) > 4:
            seq_id = parts[4]
        else:
            continue
        seq_ids.add(seq_id)

    seq_ids = sorted(seq_ids)

    Entrez.email = "barbara.r.terlouw@gmail.com"

    with open(out_file, 'w') as out:
        for chunk in chunked_iterable(seq_ids, chunk_size):
            try:
                handle = Entrez.efetch(db="protein", id=chunk, rettype="gb", retmode="text")
                for seq_record in SeqIO.parse(handle, "gb"):
                    taxonomy = seq_record.annotations.get("taxonomy", [])
                    line = f"{seq_record.id}\t" + "\t".join(taxonomy)
                    out.write(f"{line}\n")
                    print(line)
                handle.close()
                time.sleep(delay)  # Be polite to NCBI servers
            except HTTPError as e:
                print(f"[WARNING] Failed to fetch chunk {chunk}: {e}")

def get_mibig_taxonomy(mibig_fasta, parasect_fasta, out_dir):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    taxonomy = os.path.join(out_dir, "taxonomy.txt")
    fungal = os.path.join(out_dir, "fungal.fasta")
    bacterial = os.path.join(out_dir, "bacterial.fasta")
    unknown_tax = os.path.join(out_dir, "unknown_taxonomy.fasta")
    other_tax = os.path.join(out_dir, "other_taxonomy.fasta")
    in_parasect = os.path.join(out_dir, "in_parasect.fasta")
    parasect_sequences = read_fasta(parasect_fasta)

    if not os.path.exists(taxonomy):
        get_taxid(mibig_fasta, taxonomy)
    protein_to_taxonomy = parse_taxonomy(taxonomy)
    domain_to_seq = read_fasta(mibig_fasta)
    with open(fungal, 'w') as fungal_out:
        with open(bacterial, 'w') as bacterial_out:
            with open(unknown_tax, 'w') as unknown_out:
                with open(other_tax, 'w') as other_tax_out:
                    with open(in_parasect, 'w') as in_parasect_out:
                        for domain, seq in domain_to_seq.items():
                            for parasect_domain, seq_2 in parasect_sequences.items():
                                if sequences_are_equivalent(seq, seq_2):
                                    in_parasect_out.write(f">{domain}\n{seq}\n")
                                    break
                            else:
                                protein_id = domain.split('|')[4]
                                if protein_id in protein_to_taxonomy:
                                    tax = protein_to_taxonomy[protein_id]
                                    if 'Fungi' in tax:
                                        fungal_out.write(f">{domain}\n{seq}\n")
                                    elif 'Bacteria' in tax:
                                        bacterial_out.write(f">{domain}\n{seq}\n")
                                    else:
                                        other_tax_out.write(f">{domain}\n{seq}\n")
                                else:
                                    unknown_out.write(f">{domain}\n{seq}\n")


if __name__ == "__main__":
    get_mibig_taxonomy(argv[1], argv[2], argv[3])
