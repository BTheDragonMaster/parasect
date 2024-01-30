from argparse import ArgumentParser

import blosum as bl

from paras.scripts.parsers.fasta import read_fasta

BLOSUM62 = dict(bl.BLOSUM(62))


def parse_arguments():
    parser = ArgumentParser(description="Compare active sites")
    parser.add_argument("-f", required=True, type=str, help="Input fasta file containing active site signatures")
    parser.add_argument("-i", required=True, type=str,
                        help="String of active site signature to compare. Has to be of same length as signatures in fasta")
    parser.add_argument("-m", default='min', type=str,
                        help="Accepts 'max' or 'min': 'max' for most and 'min' for least similar sequence.")
    parser.add_argument("-o", required=True, type=str,
                        help="Output file")

    arguments = parser.parse_args()
    return arguments


def run():
    args = parse_arguments()

    assert args.m in ('min', 'max')
    id_to_seq = read_fasta(args.f)
    if args.m == 'min':
        score = 10000000.0
    elif args.m == 'max':
        score = -10000000.0
    else:
        raise ValueError(f"Expected 'min' or 'max'. Got '{args.m}'")

    protein_id = None

    ids_and_scores = []

    for seq_id, sequence in id_to_seq.items():
        seq_score = compare_active_sites(args.i, sequence)
        ids_and_scores.append((seq_id, seq_score))

        if args.m == 'min':
            if seq_score < score:
                score = seq_score
                protein_id = seq_id
        elif args.m == 'max':
            if seq_score > score:
                score = seq_score
                protein_id = seq_id
        else:
            raise ValueError(f"Expected 'min' or 'max'. Got '{args.m}'")

    ids_and_scores.sort(key=lambda x: x[1])
    with open(args.o, 'w') as out:
        for seq_id, seq_score in ids_and_scores:
            out.write(f"{seq_id}\t{seq_score:.1f}\n")

    return protein_id, score


def compare_active_sites(sequence_1, sequence_2):
    assert len(sequence_1) == len(sequence_2)
    score = 0.0
    for i, char_1 in enumerate(sequence_1):
        char_2 = sequence_2[i]
        score += BLOSUM62[char_1][char_2]
    return score


if __name__ == "__main__":
    print(run())
