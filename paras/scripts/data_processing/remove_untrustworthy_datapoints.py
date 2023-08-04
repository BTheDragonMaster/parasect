from sys import argv
from paras.scripts.data_processing.sequence_processing.sort_sequence_per_substrate import remove_gaps
from paras.scripts.parsers.fasta import read_fasta


def parse_untrustworthy(untrustworthy_datapoints):
    domain_ids = []

    with open(untrustworthy_datapoints, 'r') as untrustworthy:
        for line in untrustworthy:
            line = line.strip()
            domain_id, reason = line.split('\t')
            if domain_id not in domain_ids:
                domain_ids.append(domain_id)
            else:
                print(domain_id)

    return domain_ids


def remove_untrustworthy_from_alignment(alignment_fasta, untrustworthy_datapoints, out_path):
    domains_to_remove = parse_untrustworthy(untrustworthy_datapoints)
    domains = []
    sequences = []
    domain_to_seq = read_fasta(alignment_fasta)
    for domain, seq in domain_to_seq.items():
        if domain not in domains_to_remove:
            domains.append(domain)
            sequences.append(seq)

    sequences = remove_gaps(sequences)

    with open(out_path, 'w') as out:
        for i, domain in enumerate(domains):
            out.write(f">{domain}\n{sequences[i]}\n")


def remove_untrustworthy_datapoints(parasect_data, untrustworthy_datapoints, out_path):

    domain_ids = parse_untrustworthy(untrustworthy_datapoints)
    removed_domains = []

    with open(parasect_data, 'r') as parasect:
        with open(out_path, 'w') as out:
            header = parasect.readline()
            out.write(header)
            for line in parasect:
                domain_id = line.split('\t')[0]
                if domain_id not in domain_ids:
                    out.write(line)
                else:
                    removed_domains.append(domain_id)

    print(set(domain_ids).difference(set(removed_domains)))


if __name__ == "__main__":
    # remove_untrustworthy_datapoints(argv[1], argv[2], argv[3])
    remove_untrustworthy_from_alignment(argv[1], argv[2], argv[3])
