import math
from sys import argv
import time

from Bio.Align import substitution_matrices
from Bio.Align import PairwiseAligner

from paras.scripts.parsers.fasta import read_fasta
from paras.scripts.parsers.parsers import parse_domain_list

LAMBDA = 0.252
K = 0.035


def get_pairwise_identity(query: str, subject: str, aligner: PairwiseAligner):
    alignment = aligner.align(query, subject)[0]
    aln1, aln2 = alignment
    identities = len([i for i, k in zip(aln1, aln2) if i == k and i != '-'])
    identity = 100 * identities / len(aln1)
    return identity


def calculate_scores(fasta_file, out_file):
    domain_to_seq = read_fasta(fasta_file)
    domains = list(domain_to_seq.keys())
    domains.sort()

    with open(out_file, 'w') as out:

        for i, domain_1 in enumerate(domains):
            for j, domain_2 in enumerate(domains):
                if j > i:
                    aligner = PairwiseAligner()

                    aligner.mode = 'global'
                    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
                    aligner.open_gap_score = -5.0
                    aligner.extend_gap_score = -1.0
                    aligner.target_end_gap_score = 0.0
                    aligner.query_end_gap_score = 0.0

                    score = aligner.score(domain_to_seq[domain_1].upper(), domain_to_seq[domain_2].upper())
                    bitscore = (LAMBDA * score - math.log(K, math.e))/math.log(2, math.e)
                    out.write(f"{domain_1}\t{domain_2}\t{bitscore}\n")


def find_closest_identity(fasta_file, test_domain_list, train_domain_list, out_file):
    test_domains = parse_domain_list(test_domain_list)
    train_domains = parse_domain_list(train_domain_list)
    domain_to_seq = read_fasta(fasta_file)
    domains = list(domain_to_seq.keys())
    domains.sort()

    overall_highest_identity = 0.0

    with open(out_file, 'w') as out:

        for i, domain_1 in enumerate(test_domains):
            highest_identity = 0.0
            for domain_2 in train_domains:
                aligner = PairwiseAligner()

                aligner.mode = 'local'
                aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
                aligner.open_gap_score = -5.0
                aligner.extend_gap_score = -1.0
                aligner.target_end_gap_score = 0.0
                aligner.query_end_gap_score = 0.0

                identity = get_pairwise_identity(domain_to_seq[domain_1].upper(),
                                                 domain_to_seq[domain_2].upper(), aligner)
                if identity > highest_identity:
                    highest_identity = identity

                if identity > overall_highest_identity:
                    overall_highest_identity = identity

            print(f"Processed {domain_1}. Highest overall identity: {highest_identity}")

            out.write(f"{domain_1}\t{highest_identity:.2f}\n")

    print(overall_highest_identity)


def pairwise_identity(fasta_file, domain_list, out_file):
    domains = parse_domain_list(domain_list)
    domain_to_seq = read_fasta(fasta_file)

    with open(out_file, 'w') as out:

        for i, domain_1 in enumerate(domains):
            start_time = time.time()
            for j, domain_2 in enumerate(domains):
                if j > i:
                    aligner = PairwiseAligner()

                    aligner.mode = 'local'
                    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
                    aligner.open_gap_score = -5.0
                    aligner.extend_gap_score = -1.0
                    aligner.target_end_gap_score = 0.0
                    aligner.query_end_gap_score = 0.0

                    identity = get_pairwise_identity(domain_to_seq[domain_1].upper(),
                                                     domain_to_seq[domain_2].upper(), aligner)



                    out.write(f"{domain_1}\t{domain_2}\t{identity}\n")
            end_time = time.time()
            print(f"Processed {domain_1}. Time: {end_time - start_time}")


if __name__ == "__main__":
    # calculate_scores(argv[1], argv[2])
    # find_closest_identity(argv[1], argv[2], argv[3], argv[4])
    pairwise_identity(argv[1], argv[2], argv[3])

