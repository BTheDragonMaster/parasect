from sys import argv
from paras.scripts.parsers.fasta import read_fasta


def compare_active_site_files(file_1, file_2):
    id_to_seq_1 = read_fasta(file_1)
    id_to_seq_2 = read_fasta(file_2)
    matching = 0
    mismatching = 0
    gaps_1 = 0
    gaps_2 = 0
    for id_1, seq_1 in id_to_seq_1.items():
        if id_1 in id_to_seq_2:
            seq_2 = id_to_seq_2[id_1]
            gaps_1 += seq_1.count('-')
            gaps_2 += seq_2.count('-')
            # if '-' in seq_1 or '-' in seq_2:
            #     print(f"{id_1}\t{seq_1}\t{seq_2}")
            if seq_1 != seq_2:
                print(f"{id_1}\t{seq_1}\t{seq_2}")
                mismatching += 1
                # if '-' in seq_2:
                #     print(f"{id_1}\t{seq_1}\t{seq_2}")
            else:
                matching += 1
        else:
            print(f"{id_1} not found")

    print(f"Matching sequences: {matching}")
    print(f"Mismatching sequences: {mismatching}")
    print(f"Gaps file 1: {gaps_1}")
    print(f"Gaps file 2: {gaps_2}")


if __name__ == "__main__":
    compare_active_site_files(argv[1], argv[2])
