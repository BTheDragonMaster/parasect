import os
from sys import argv

from paras.scripts.feature_extraction.sequence_feature_extraction.hmm.run_hmmscan import run_hmmscan, HMM3_FILE
from paras.scripts.feature_extraction.sequence_feature_extraction.extract_domains import parse_hmm3_results, \
    hits_to_domains


def get_domains(fasta_file, out_dir):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    hmm_file = os.path.join(out_dir, "hmm_output.hmm")
    active_site_10_out_file = os.path.join(out_dir, "active_site_10.fasta")
    active_site_34_out_file = os.path.join(out_dir, "active_site_34.fasta")
    sequences_out_file = os.path.join(out_dir, "domain_sequences.fasta")
    run_hmmscan(HMM3_FILE, fasta_file, hmm_file)
    id_to_hit = parse_hmm3_results(hmm_file)
    domains = hits_to_domains(id_to_hit, fasta_file, profile=True)
    with open(active_site_10_out_file, 'w') as active_site_10_out:
        with open(active_site_34_out_file, 'w') as active_site_34_out:
            with open(sequences_out_file, 'w') as sequences_out:
                for domain in domains:
                    domain.domain_id = f"{domain.protein_name}.A{domain.domain_nr}"
                    domain.write_sequence(active_site_10_out, "signature")
                    domain.write_sequence(active_site_34_out, "extended_signature")
                    domain.write_sequence(sequences_out, "full")


if __name__ == "__main__":
    get_domains(argv[1], argv[2])
