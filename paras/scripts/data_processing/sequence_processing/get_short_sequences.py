from sys import argv

from paras.scripts.parsers.fasta import read_fasta, write_fasta

limit = int(argv[3])

id_to_seq = read_fasta(argv[1])
short_id_to_seq = {}
for seq_id, seq in id_to_seq.items():
    if len(seq) <= limit:
        short_id_to_seq[seq_id] = seq

write_fasta(short_id_to_seq, argv[2])
