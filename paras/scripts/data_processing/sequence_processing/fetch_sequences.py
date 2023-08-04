import requests
from ratelimiter import RateLimiter
from io import StringIO

from paras.scripts.parsers.fasta import read_fasta
from paras.scripts.parsers.tabular import Tabular

MAX_CALLS = 3


@ RateLimiter(MAX_CALLS)
def get_from_ncbi(accession: str) -> str:
    params = {
        "db": "protein",
        "rettype": "fasta",
        "id": accession,
    }
    res = requests.get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", params)

    if res.status_code != 200:
        print(f"Could not find sequence for {accession}.")
        return ''

    return res.text


def fetch_sequences(specificities_file, out_file):
    protein_ids = set()
    mibig_data = Tabular(specificities_file, [2, 4])
    for seq_id in mibig_data.data:
        protein_id = mibig_data.data[seq_id]["protein_name"]
        protein_ids.add(protein_id)

    id_to_seq = {}

    for protein_id in protein_ids:
        if protein_id:
            fasta = StringIO(get_from_ncbi(protein_id))
            if fasta:
                try:
                    id_to_sequence = read_fasta(fasta)
                    sequence = None
                    for seq_id, seq in id_to_sequence.items():
                        sequence = seq
                        break
                    if sequence:
                        id_to_seq[protein_id] = sequence
                    else:
                        print(f"Could not find sequence for {protein_id}.")
                except Exception as e:
                    print(e)
                    continue

    with open(out_file, 'w') as out:
        for prot_id, seq in id_to_seq.items():
            out.write(f">{prot_id}\n{seq}\n")
