import os
import paras.data.sequence_data.hmm
import paras.data.sequence_data.sequences

from paras.scripts.feature_extraction.sequence_feature_extraction.hmm.run_hmmscan import run_hmmscan

REF_SEQ_FILE = os.path.join(os.path.dirname(paras.data.sequence_data.sequences.__file__), 'reference_sequence.fasta')
HMM_FILE_SEQUENCE = os.path.join(os.path.dirname(paras.data.sequence_data.hmm.__file__), 'AMP-binding_full.hmm')
HMM_FILE_HYBRID = os.path.join(os.path.dirname(paras.data.sequence_data.hmm.__file__), 'AMP-binding_hybrid.hmm')
HMM_FILE_STRUCTURE = os.path.join(os.path.dirname(paras.data.sequence_data.hmm.__file__), 'AMP-binding_structure.hmm')

if __name__ == "__main__":
    run_hmmscan(HMM_FILE_SEQUENCE, REF_SEQ_FILE, 'sequence.hmm_result')
    run_hmmscan(HMM_FILE_STRUCTURE, REF_SEQ_FILE, 'structure.hmm_result')
    run_hmmscan(HMM_FILE_HYBRID, REF_SEQ_FILE, 'hybrid.hmm_result')

