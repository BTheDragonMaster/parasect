import os
from sys import argv

from Bio import Entrez
from Bio import SeqIO

from paras.scripts.parsers.fasta import read_fasta
import paras.data.sequence_data.sequences

DOMAINS = os.path.join(os.path.dirname(paras.data.sequence_data.sequences.__file__), "domains.fasta")


def retrieve_gbk_record(fasta_file, out_file):
    id_to_seq = read_fasta(fasta_file)
    seq_ids = []

    for seq_id_str in id_to_seq.keys():
        seq_id_list = seq_id_str.split('|')
        for seq_id in seq_id_list:
            seq_ids.append('.'.join(seq_id.split('.')[:-1]))

    seq_ids = list(set(seq_ids))
    seq_ids.sort()

    Entrez.email = "barbara.r.terlouw@gmail.com"
    handle = Entrez.efetch(db="protein", id=seq_ids, rettype="gb", retmode="text")
    with open(out_file, 'w') as out:
        for seq_record in SeqIO.parse(handle, "gb"):
            line = f"{seq_record.id}\t" + "\t".join(seq_record.annotations["taxonomy"])
            out.write(f"{line}\n")
            print(line)
    handle.close()


if __name__ == "__main__":
    retrieve_gbk_record(DOMAINS, argv[1])