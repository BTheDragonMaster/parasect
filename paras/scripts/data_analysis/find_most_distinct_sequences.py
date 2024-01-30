from argparse import ArgumentParser

import blosum as bl

from paras.scripts.parsers.fasta import read_fasta
from paras.scripts.parsers.parsers import parse_taxonomy

BLOSUM62 = dict(bl.BLOSUM(62))


def parse_arguments():
    parser = ArgumentParser(description="Compare active sites")
    parser.add_argument("-f", required=True, type=str, help="Input fasta file containing active site signatures")
    parser.add_argument("-t", default=None, type=str, help="Input taxonomy file")
    parser.add_argument("-l", default=None, type=float, help="Threshold below which to report dissimilar sequences")
    arguments = parser.parse_args()
    return arguments


def run():
    args = parse_arguments()


    id_to_seq = read_fasta(args.f)
    min_score = 10000000.0
    max_score = -100000000.0
    most_distinct_pair = None
    most_similar_pair = None

    seq_ids = list(id_to_seq.keys())
    seq_ids.sort()

    print(len(seq_ids))

    if args.t is not None:
        id_to_taxonomy = parse_taxonomy(args.t)
        for seq_id in seq_ids[:]:
            if seq_id not in id_to_taxonomy:
                seq_ids.remove(seq_id)
            elif id_to_taxonomy[seq_id][0] != 'Bacteria':
                seq_ids.remove(seq_id)

    print(len(seq_ids))

    dissimilar_pairs = []

    for i, seq_id_1 in enumerate(seq_ids):
        seq_1 = id_to_seq[seq_id_1]
        for seq_id_2 in seq_ids[i + 1:]:
            seq_2 = id_to_seq[seq_id_2]
            if '-' not in seq_1 and '-' not in seq_2:
                seq_score = compare_active_sites(seq_1, seq_2)
                if args.l is not None and seq_score < args.l:
                    dissimilar_pairs.append((seq_id_1, seq_id_2, seq_score))
                if seq_score < min_score:
                    min_score = seq_score
                    most_distinct_pair = (seq_id_1, seq_id_2)
                if seq_score > max_score:
                    max_score = seq_score
                    most_similar_pair = (seq_id_1, seq_id_2)

    print("Most similar pair:")
    print(f"{most_similar_pair[0]}\t{most_similar_pair[1]}\t{max_score:.1f}")

    print("Most distinct pair:")
    print(f"{most_distinct_pair[0]}\t{most_distinct_pair[1]}\t{min_score:.1f}")

    dissimilar_pairs.sort(key=lambda x: x[2])
    print("Dissimilar pairs:")
    for id_1, id_2, seq_score in dissimilar_pairs:
        print(f"{id_1}\t{id_2}\t{seq_score:.2f}")


def compare_active_sites(sequence_1, sequence_2):
    assert len(sequence_1) == len(sequence_2)
    score = 0.0
    for i, char_1 in enumerate(sequence_1):

        char_2 = sequence_2[i]
        score += BLOSUM62[char_1][char_2]
    return score


if __name__ == "__main__":
    run()
