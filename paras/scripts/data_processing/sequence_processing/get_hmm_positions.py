#!/usr/bin/env python

import os

from paras.scripts.feature_extraction.sequence_feature_extraction.read_positions import POSITIONS_EXTENDED_SIGNATURE, \
    POSITIONS_SIGNATURE
from paras.scripts.feature_extraction.sequence_feature_extraction.extract_domains import parse_hmm_results, \
    parse_domain_id, parse_hmm2_results
from paras.scripts.feature_extraction.sequence_feature_extraction.hmm.run_hmmscan import run_hmmscan, run_hmmpfam2
import paras.data.sequence_data.hmm
import paras.data.sequence_data.sequences

REF_SEQ_FILE = os.path.join(os.path.dirname(paras.data.sequence_data.sequences.__file__), 'reference_sequence.fasta')
HMM_FILE_SEQUENCE = os.path.join(os.path.dirname(paras.data.sequence_data.hmm.__file__), 'AMP-binding_full.hmm')
HMM_FILE_HYBRID = os.path.join(os.path.dirname(paras.data.sequence_data.hmm.__file__), 'AMP-binding_hybrid.hmm')
HMM_FILE_STRUCTURE = os.path.join(os.path.dirname(paras.data.sequence_data.hmm.__file__), 'AMP-binding_structure.hmm')
HMM_FILE_NRPSPREDICTOR = os.path.join(os.path.dirname(paras.data.sequence_data.hmm.__file__), 'AMP-binding_nrpspredictor.hmm')
HMM_FILE_HMMER2 = os.path.join(os.path.dirname(paras.data.sequence_data.hmm.__file__), 'AMP-binding_hmmer2.hmm')


def find_hmm_positions(id_to_hit, original_positions):
    for seq_id, hit in id_to_hit.items():
        sequence_id, hit_id, hit_start, hit_end = parse_domain_id(seq_id)
        if hit_id == 'AMP-binding':
            print(hit.aln[0].seq)
            print(hit.aln[1].seq)
            hmm_seq = hit.aln[1].seq
            query_seq = hit.aln[0].seq
            hmm_start = hit.hit_start
            seq_start = hit.query_start

            positions_hmm = []
            position_hmm = hmm_start
            position_query = seq_start
            for i, amino in enumerate(hmm_seq):
                query_amino = query_seq[i]
                if 160 < i < 200:
                    print(i, position_hmm, position_query, query_amino, amino)

                if query_amino != '-':
                    if position_query in original_positions:
                        positions_hmm.append(position_hmm)

                    position_query += 1

                if amino not in "-.":
                    position_hmm += 1

            return positions_hmm

    return None


if __name__ == "__main__":
    # run_hmmscan(HMM_FILE_SEQUENCE, REF_SEQ_FILE, 'sequence.hmm_result')
    # run_hmmscan(HMM_FILE_STRUCTURE, REF_SEQ_FILE, 'structure.hmm_result')
    # run_hmmscan(HMM_FILE_HYBRID, REF_SEQ_FILE, 'hybrid.hmm_result')
    run_hmmscan(HMM_FILE_NRPSPREDICTOR, REF_SEQ_FILE, 'nrpspredictor.hmm_result')
    run_hmmpfam2(HMM_FILE_HMMER2, REF_SEQ_FILE, 'hmmer2.hmm_result')

    sequence_hit_to_id = parse_hmm_results('sequence.hmm_result')
    nrpspredictor_hit_to_id = parse_hmm_results('nrpspredictor.hmm_result')
    hmmer_2_hit_to_id = parse_hmm2_results('hmmer2.hmm_result')
    # structure_hit_to_id = parse_hmm_results('structure.hmm_result')
    # hybrid_hit_to_id = parse_hmm_results('hybrid.hmm_result')

    # print('\t'.join(map(str, find_hmm_positions(sequence_hit_to_id, POSITIONS_EXTENDED_SIGNATURE))))
    # print('\t'.join(map(str, find_hmm_positions(sequence_hit_to_id, POSITIONS_SIGNATURE))))

    print('\nHmmer3:')
    print('\t'.join(map(str, find_hmm_positions(nrpspredictor_hit_to_id, POSITIONS_EXTENDED_SIGNATURE))))
    print('\t'.join(map(str, find_hmm_positions(nrpspredictor_hit_to_id, POSITIONS_SIGNATURE))))

    print('\nHmmer2:')
    print('\t'.join(map(str, find_hmm_positions(hmmer_2_hit_to_id, POSITIONS_EXTENDED_SIGNATURE))))
    print('\t'.join(map(str, find_hmm_positions(hmmer_2_hit_to_id, POSITIONS_SIGNATURE))))

