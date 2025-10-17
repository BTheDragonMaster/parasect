from Bio import Entrez
import json
import os
import time
from urllib.error import HTTPError
from sys import argv

from paras.scripts.parsers.fasta import read_fasta
from paras.scripts.parsers.parsers import parse_taxonomy


def chunked_iterable(iterable, size):
    for i in range(0, len(iterable), size):
        yield iterable[i:i + size]


def get_lineage_names(record):
    """
    Extract a list of lineage names from a single taxonomy record.
    """
    lineage = [item['ScientificName'] for item in record['LineageEx']]
    lineage.append(record['ScientificName'])  # Add the organism itself at the end
    return lineage


def get_taxid(fasta_file, mibig_dir, out_file, chunk_size=200, delay=0.4):
    id_to_seq = read_fasta(fasta_file)

    id_to_ncbi_tax = {}
    tax_ids = set()

    # Extract tax IDs from the JSON metadata
    for seq_id_str in id_to_seq.keys():
        parts = seq_id_str.split('|')
        if len(parts) > 0:
            mibig_id = parts[0].split('.')[0]
        else:
            continue
        json_file = os.path.join(mibig_dir, f"{mibig_id}.json")
        with open(json_file, 'r') as json_data:
            cluster_data = json.load(json_data)
            try:
                tax_id = cluster_data["taxonomy"]["ncbiTaxId"]
                id_to_ncbi_tax[seq_id_str] = tax_id
                tax_ids.add(tax_id)
            except KeyError:
                print(f"Couldn't find tax_id for {mibig_id}")
                continue

    tax_ids = sorted(tax_ids)

    Entrez.email = "barbara.r.terlouw@gmail.com"

    taxid_to_lineage = {}


    for chunk in chunked_iterable(tax_ids, chunk_size):
        try:
            handle = Entrez.efetch(db="taxonomy", id=chunk, retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            time.sleep(delay)

            for record in records:
                taxid = int(record['TaxId'])
                lineage_names = get_lineage_names(record)
                taxid_to_lineage[taxid] = lineage_names

        except HTTPError as e:
            print(f"[WARNING] Failed to fetch chunk {chunk}: {e}")

    with open(out_file, 'w') as out:
        for seq_id, tax_id in id_to_ncbi_tax.items():
            if tax_id in taxid_to_lineage:
                taxonomy = taxid_to_lineage[tax_id]
                line = f"{seq_id}\t" + "\t".join(taxonomy)
                out.write(f"{line}\n")
            else:
                print(f"Couldn't find taxonomy for {seq_id}")


def extend_existing(tax_file, unknown_file, bacterial_file, fungal_file, other_file):
    taxonomy = parse_taxonomy(tax_file)
    unknowns = read_fasta(unknown_file)
    with open(bacterial_file, 'a') as bacterial:
        with open(fungal_file, 'a') as fungal:
            with open(other_file, 'a') as other:
                for seq_id, seq in unknowns.items():
                    tax = taxonomy[seq_id]
                    if 'Fungi' in tax:
                        fungal.write(f">{seq_id}\n{seq}\n")
                    elif 'Bacteria' in tax:
                        bacterial.write(f">{seq_id}\n{seq}\n")
                    else:
                        other.write(f">{seq_id}\n{seq}\n")


if __name__ == "__main__":
    extend_existing(argv[1], argv[2], argv[3], argv[4], argv[5])
