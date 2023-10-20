import subprocess
import os
import paras.data.sequence_data.hmm


HMM_FILE = os.path.join(os.path.dirname(paras.data.sequence_data.hmm.__file__), 'AMP-binding_full.hmm')


def run_hmmscan(hmm_dir, fasta_file, hmm_out):
    """
    Run hmmscan from command line

    Input:
    hmm_dir: str, dir of .hmm file containing the HMMs to be used in the scan
    fasta_dir: str, dir of .fasta file containing the sequences to be scanned
    out_dir: str, file location containing results of hmmscan

    """

    with open(hmm_out, 'w') as out:
        command = ['hmmscan', hmm_dir, fasta_file]
        subprocess.call(command, stdout=out)
