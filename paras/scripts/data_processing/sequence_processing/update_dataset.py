from compare_domains import sequences_are_equivalent
from paras.scripts.parsers.tabular import Tabular
from paras.scripts.parsers.fasta import read_fasta

from sys import argv


def update_sequences(domains_in, domains_new, dataset_in, dataset_out):
    id_to_seq_old = read_fasta(domains_in)
    id_to_seq_new = read_fasta(domains_new)
    parasect_data = Tabular(dataset_in, [0])
    data_to_seq = {}

    all_ids = set()
    id_to_full = {}

    for seq_id_1, seq_1 in id_to_seq_old.items():
        if seq_id_1 not in parasect_data.get_column("domain_id"):
            print("Not in dataset", seq_id_1)
        seq_ids = seq_id_1.split('|')
        for seq_id in seq_ids:
            id_to_full[seq_id] = seq_id_1
            if seq_id in all_ids:
                print(f"Double id: {seq_id}")
            else:
                all_ids.add(seq_id)

    missing_ids = []

    for seq_id_1, seq_1 in id_to_seq_old.items():
        seq_ids = seq_id_1.split('|')
        match_found = False
        for seq_id in seq_ids:
            seq_id_found = False
            for seq_id_2, seq_2 in id_to_seq_new.items():

                if seq_id == seq_id_2:
                    seq_id_found = True
                    match_found = True
                    assert sequences_are_equivalent(seq_1.upper(), seq_2.upper())
                    new_seq = seq_2.upper()
                    data_to_seq[seq_id_1] = new_seq
                    break

            if not seq_id_found:
                missing_ids.append(seq_id)

        if not match_found:
            raise ValueError(f"Couldn't find sequence for {seq_id_1}")

    renamed_datapoints = {}

    for missing_id in missing_ids:
        old_full_name = id_to_full[missing_id]
        new_full_name = old_full_name.split('|')
        new_full_name.remove(missing_id)
        new_full_name = '|'.join(new_full_name)
        renamed_datapoints[old_full_name] = new_full_name

    with open(dataset_out, 'w') as out:
        out.write(f"domain_id\tsequence\tspecificity\n")
        for datapoint in parasect_data.data:
            domain_id = parasect_data.get_value(datapoint, "domain_id")
            sequence = data_to_seq[domain_id]
            specificity = parasect_data.get_value(datapoint, "specificity")
            if domain_id in renamed_datapoints:
                domain_id = renamed_datapoints[domain_id]

            out.write(f"{domain_id}\t{sequence}\t{specificity}\n")

    clean_data = Tabular(dataset_out, [0])
    for datapoint in clean_data.data:
        sequence = clean_data.get_value(datapoint, "sequence")
        if len(sequence) < 300:
            print(clean_data.get_value(datapoint, "domain_id"), len(sequence))

    sequences = clean_data.get_column("sequence")
    print(min([len(s) for s in sequences]), max([len(s) for s in sequences]))

    print(len([seq for seq in sequences if len(seq) < 400]))


if __name__ == "__main__":
    update_sequences(argv[1], argv[2], argv[3], argv[4])