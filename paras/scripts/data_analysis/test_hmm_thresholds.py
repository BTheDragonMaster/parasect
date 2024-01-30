from sys import argv
import os

from Bio.SearchIO._model import HSP
from Bio import SearchIO

from paras.scripts.feature_extraction.sequence_feature_extraction.sequence_labels import make_domain_id
from paras.scripts.feature_extraction.sequence_feature_extraction.extract_domains import hits_to_domains
from paras.scripts.parsers.fasta import read_fasta
from paras.scripts.feature_extraction.sequence_feature_extraction.hmm.run_hmmscan import HMM2_FILE, \
    run_hmmpfam2
from paras.scripts.data_processing.temp import TEMP_DIR, clear_temp
from paras.scripts.feature_extraction.sequence_feature_extraction.rename_sequences import rename_sequences, parse_mapping


def test_hmm2_thresholds(fasta_in, hmm_out, signature_out_folder):
    if not os.path.exists(signature_out_folder):
        os.mkdir(signature_out_folder)
    stach_out = os.path.join(signature_out_folder, "active_site_10.fasta")
    sig_out = os.path.join(signature_out_folder, "active_site_34.fasta")
    mapping, renamed_fasta_file = rename_sequences(fasta_in, TEMP_DIR)
    new_to_original = parse_mapping(mapping)
    all_ids = read_fasta(renamed_fasta_file).keys()
    nr_domains = len(all_ids)
    run_hmmpfam2(HMM2_FILE, renamed_fasta_file, os.path.join(TEMP_DIR, hmm_out))
    threshold_to_id_to_hit: dict[float, dict[str, HSP]] = {}

    all_ids = set(all_ids)

    for result in SearchIO.parse(os.path.join(TEMP_DIR, hmm_out), 'hmmer2-text'):
        for hsp in result.hsps:
            for threshold in [0.001]:
                if threshold not in threshold_to_id_to_hit:
                    threshold_to_id_to_hit[threshold]: dict[str, HSP] = {}
                if hsp.evalue < threshold:
                    if hsp.hit_id == 'AMP-binding' or hsp.hit_id == 'AMP-binding_C':
                        header = make_domain_id(result.id, hsp.hit_id, hsp.query_start, hsp.query_end)
                        threshold_to_id_to_hit[threshold][header] = hsp

    print(nr_domains)
    sigs = set()

    with open(stach_out, 'w') as stach:
        with open(sig_out, 'w') as sig:
            for threshold, id_to_hit in threshold_to_id_to_hit.items():
                a_domains = hits_to_domains(id_to_hit, renamed_fasta_file)
                for a_domain in a_domains:
                    stach.write(f">{new_to_original[a_domain.protein_name]}\n{a_domain.signature}\n")
                    sig.write(f">{new_to_original[a_domain.protein_name]}\n{a_domain.extended_signature}\n")
                    sigs.add(a_domain.extended_signature)
                    if a_domain.protein_name in all_ids:
                        all_ids.remove(a_domain.protein_name)

                print(threshold, len(a_domains))

    print(len(sigs))

    for seq_id in all_ids:
        print(new_to_original[seq_id])



    clear_temp()


if __name__ == "__main__":
    test_hmm2_thresholds(argv[1], argv[2], argv[3])

