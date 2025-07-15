import os
from sys import argv

from paras.scripts.feature_extraction.sequence_feature_extraction.hmm.run_hmmscan import run_hmmscan, HMM3_FILE
from paras.scripts.feature_extraction.sequence_feature_extraction.extract_domains import parse_hmm3_results
from paras.scripts.feature_extraction.sequence_feature_extraction.rename_sequences import rename_sequences, parse_mapping
from paras.scripts.parsers.fasta import read_fasta, write_fasta
from paras.scripts.feature_extraction.sequence_feature_extraction.sequence_labels import parse_domain_id
from paras.scripts.feature_extraction.sequence_feature_extraction.adenylation_domain import AdenylationDomain


def hits_to_domains(id_to_hit, fasta_file, out_file, mapping_file):
    hits_by_seq_id = {}

    for hit_key, hit in id_to_hit.items():
        seq_id, hit_id, hit_start, hit_end = parse_domain_id(hit_key)
        if seq_id not in hits_by_seq_id:
            hits_by_seq_id[seq_id] = []

        hits_by_seq_id[seq_id].append((hit_id, hit_start, hit_end, hit_key))

    counter = 0

    seq_id_to_domains = {}

    seq_id_with_match = []
    seq_id_without_match = []
    seq_id_without_hits = []
    seq_id_to_hit_nr = {}
    seq_id_has_c = []
    seq_id_has_n = []

    for seq_id, hits in hits_by_seq_id.items():

        counter += 1
        seq_id_to_hit_nr[seq_id] = len(hits)
        if not hits:
            seq_id_without_hits.append(seq_id)
        for hit_id_1, hit_start_1, hit_end_1, hit_key_1 in hits:
            if hit_id_1 == 'AMP-binding':
                seq_id_has_n.append(seq_id)
                if seq_id not in seq_id_to_domains:
                    seq_id_to_domains[seq_id] = []
                match_found = False
                for hit_id_2, hit_start_2, hit_end_2, hit_key_2 in hits:
                    if hit_id_2 == 'AMP-binding_C':
                        seq_id_has_c.append(seq_id)
                        # Todo: check that 200 is a good cutoff score

                        if (hit_start_2 > hit_end_1) and (hit_start_2 - hit_end_1 < 200):
                            a_domain = AdenylationDomain(seq_id, hit_start_1, hit_end_2)
                            seq_id_to_domains[seq_id].append(a_domain)
                            match_found = True
                            seq_id_with_match.append(seq_id)
                            break

                if not match_found:
                    seq_id_without_match.append(seq_id)
                    a_domain = AdenylationDomain(seq_id, hit_start_1, hit_end_1)
                    a_domain.set_domain_signatures_hmm(id_to_hit[hit_key_1])
                    seq_id_to_domains[seq_id].append(a_domain)

    for seq_id, domains in seq_id_to_domains.items():
        domains.sort(key=lambda x: x.start)

    new_to_old = parse_mapping(mapping_file)

    number_without_match = 0
    number_without_match_without_c = 0
    number_without_match_without_c_long = 0

    with open(out_file, 'w') as hit_tracker:
        hit_tracker.write(f"domain_id\tnr_hits\thas_N\thas_C\thas_matching_pair\tdomain_positions\n")

        for seq_id in seq_id_to_hit_nr:
            has_match = 0
            has_n = 0
            has_c = 0
            if seq_id in seq_id_with_match:
                has_match = 1
            if seq_id in seq_id_has_c:
                has_c = 1
            if seq_id in seq_id_has_n:
                has_n = 1
            range_string = ""
            domain_length = 0

            for domain in seq_id_to_domains[seq_id]:
                range_string += f"|{domain.start}-{domain.end}|"
                if domain.end - domain.start > domain_length:
                    domain_length = domain.end - domain.start

            if not has_match:
                number_without_match += 1
                if not has_c:
                    number_without_match_without_c += 1
                    if domain_length > 300:
                        number_without_match_without_c_long += 1
            hit_tracker.write(f"{new_to_old[seq_id]}\t{seq_id_to_hit_nr[seq_id]}\t{has_n}\t{has_c}\t{has_match}\t{range_string}\n")

    fasta = read_fasta(fasta_file)

    print(number_without_match_without_c_long)

    for seq_id, sequence in fasta.items():
        counter = 1
        if seq_id in seq_id_to_domains:
            for a_domain in seq_id_to_domains[seq_id]:
                assert seq_id == a_domain.protein_name
                a_domain_sequence = sequence[a_domain.start:a_domain.end]
                if len(a_domain_sequence) > 100:
                    a_domain.set_sequence(a_domain_sequence)
                    a_domain.set_domain_number(counter)
                    counter += 1

    filtered_a_domains = []

    for seq_id, a_domains in seq_id_to_domains.items():
        for a_domain in a_domains:
            if a_domain.sequence and a_domain.domain_nr:
                filtered_a_domains.append(a_domain)

    filtered_a_domains.sort(key=lambda x: (x.protein_name, x.start))

    return filtered_a_domains


def get_domains(fasta_file, out_dir, mapping_file):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    hmm_file = os.path.join(out_dir, "hmm_output.hmm")
    hit_tracker = os.path.join(out_dir, 'hmm_hit_summary.txt')
    sequences_out_file = os.path.join(out_dir, "domain_sequences.fasta")

    run_hmmscan(HMM3_FILE, fasta_file, hmm_file)
    id_to_hit = parse_hmm3_results(hmm_file)
    domains = hits_to_domains(id_to_hit, fasta_file, hit_tracker, mapping_file)
    protein_names = []
    for domain in domains:
        if domain.domain_nr == 2:
            print(domain.start, domain.end, domain.protein_name)
            protein_names.append(domain.protein_name)
    with open(sequences_out_file, 'w') as sequences_out:
        for domain in domains:
            if domain.protein_name in protein_names and domain.domain_nr == 1:
                print(domain.start, domain.end, domain.protein_name)
            domain.domain_id = f"{domain.protein_name}.A{domain.domain_nr}"
            domain.write_sequence(sequences_out, "full")


if __name__ == "__main__":
    mapping_file, renamed_fasta = rename_sequences(argv[1], argv[2])
    new_to_original: dict[str, str] = parse_mapping(mapping_file)

    get_domains(renamed_fasta, argv[2], mapping_file)
    full_dir = os.path.join(argv[2], "domain_sequences.fasta")
    id_to_seq = read_fasta(full_dir)
    old_to_seq = {}
    seen_olds = []
    for seq_id, seq in id_to_seq.items():
        old_id = new_to_original[seq_id.split('.')[0]]
        if seq_id.split('.')[1] == 'A2':
            print(old_id)
        else:
            old_to_seq[old_id] = seq
            seen_olds.append(old_id)

    for seq_id in read_fasta(argv[1]):
        if seq_id not in seen_olds:
            print("Not seen:", seq_id)

    write_fasta(old_to_seq, os.path.join(argv[2], "domain_sequences_fixed_ids.fasta"))
