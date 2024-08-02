from sys import argv
import os
from paras.scripts.data_processing.sequence_processing.sort_sequence_per_substrate import remove_gaps
from paras.scripts.parsers.fasta import read_fasta
from paras.scripts.parsers.parsers import parse_domain_list


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


def remove_untrustworthy_from_training_data(train_test_split_dir, untrustworthy_data_dir, out_path):
    if not os.path.exists(out_path):
        os.mkdir(out_path)

    untrustworthy_ids = parse_untrustworthy(untrustworthy_data_dir)

    for file_name in os.listdir(train_test_split_dir):
        in_file = os.path.join(train_test_split_dir, file_name)
        if file_name in ['train.txt', 'test.txt', 'test_filtered.txt'] and os.path.isfile(in_file):
            out_file = os.path.join(out_path, file_name)

            domain_ids = parse_domain_list(in_file)
            with open(out_file, 'w') as out:
                for domain_id in domain_ids:
                    if domain_id not in untrustworthy_ids:
                        out.write(f"{domain_id}\n")
        elif file_name == 'crossvalidation' and os.path.isdir(in_file):
            new_crossval_dir = os.path.join(out_path, file_name)
            if not os.path.exists(new_crossval_dir):
                os.mkdir(new_crossval_dir)

            for crossval_file_name in os.listdir(in_file):
                crossval_file = os.path.join(in_file, crossval_file_name)
                if 'crossval' in crossval_file_name and os.path.isfile(crossval_file):
                    out_file = os.path.join(new_crossval_dir, crossval_file_name)

                    domain_ids = parse_domain_list(crossval_file)
                    with open(out_file, 'w') as crossval_out:
                        for domain_id in domain_ids:
                            if domain_id not in untrustworthy_ids:
                                crossval_out.write(f"{domain_id}\n")


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

    print(set(removed_domains))


if __name__ == "__main__":
    # remove_untrustworthy_datapoints(argv[1], argv[2], argv[3])
    # remove_untrustworthy_from_alignment(argv[1], argv[2], argv[3])
    remove_untrustworthy_from_training_data(argv[1], argv[2], argv[3])
