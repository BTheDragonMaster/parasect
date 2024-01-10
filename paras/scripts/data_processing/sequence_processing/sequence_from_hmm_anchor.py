from sys import argv

from paras.scripts.parsers.fasta import read_fasta
from paras.scripts.feature_extraction.sequence_feature_extraction.hmm.run_hmmscan import run_hmmscan, HMM3_FILE
from paras.scripts.feature_extraction.sequence_feature_extraction.extract_domains import parse_hmm3_results, parse_domain_id
from paras.scripts.feature_extraction.sequence_feature_extraction.rename_sequences import rename_sequences, parse_mapping


def get_sequence_from_hmm_anchor(hmm_file, sequence_file, out_file, hmm_position=167, leading=5, trailing=27):
    mapping_file, renamed_fasta = rename_sequences(sequence_file, "renamed.fasta")

    id_to_seq = read_fasta(renamed_fasta)
    run_hmmscan(hmm_file, renamed_fasta, 'hmmer2.hmm_result')
    id_to_hit = parse_hmm3_results('hmmer2.hmm_result')

    id_to_index = {}

    for seq_id, hit in id_to_hit.items():
        sequence_id, hit_id, hit_start, hit_end = parse_domain_id(seq_id)
        if hit_id == 'AMP-binding':
            hmm_seq = hit.aln[1].seq
            query_seq = hit.aln[0].seq
            hmm_start = hit.hit_start
            seq_start = hit.query_start

            position_hmm = hmm_start
            position_query = seq_start
            for i, amino in enumerate(hmm_seq):
                query_amino = query_seq[i]

                if amino not in "-.":
                    if position_hmm == hmm_position:
                        assert query_amino == id_to_seq[sequence_id][position_query]
                        id_to_index[sequence_id] = position_query
                    position_hmm += 1

                if query_amino != '-':
                    position_query += 1

    mapping = parse_mapping(mapping_file)

    with open(out_file, 'w') as out:
        for seq_id, seq_index in id_to_index.items():
            sequence = id_to_seq[seq_id]
            extracted_sequence = sequence[max([0, seq_index - leading]) : seq_index + trailing + 1]
            original_id = mapping[seq_id]
            out.write(f">{original_id}\n{extracted_sequence}\n")


if __name__ == "__main__":
    get_sequence_from_hmm_anchor(HMM3_FILE, argv[1], argv[2])



