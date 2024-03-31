import subprocess
import os
import paras.data.sequence_data.hmm
from io import StringIO


HMM3_FILE = os.path.join(os.path.dirname(paras.data.sequence_data.hmm.__file__), 'AMP-binding_full.hmm')
HMM2_FILE = os.path.join(os.path.dirname(paras.data.sequence_data.hmm.__file__), 'AMP-binding_hmmer2.hmm')


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


def run_hmmpfam2(hmm_dir, fasta_file, hmm_out):
    """
    Run hmmscan from command line

    Input:
    hmm_dir: str, dir of .hmm file containing the HMMs to be used in the scan
    fasta_dir: str, dir of .fasta file containing the sequences to be scanned
    out_dir: str, file location containing results of hmmscan

    """

    with open(hmm_out, 'w') as out:
        command = ['hmmpfam2', hmm_dir, fasta_file]
        # command = [hmmer_path, hmm_dir, fasta_file]
        subprocess.call(command, stdout=out)


def run_hmmpfam2_single(hmm_dir, fasta_sequence, max_eval=0.1):
    """
    Run hmmscan from command line

    Input:
    hmm_dir: str, dir of .hmm file containing the HMMs to be used in the scan
    fasta_dir: str, dir of .fasta file containing the sequences to be scanned
    out_dir: str, file location containing results of hmmscan

    """

    if fasta_sequence is not None:
        stdin = subprocess.PIPE
        input_bytes = fasta_sequence.encode("utf-8")
    else:
        stdin = None
        input_bytes = None

    with subprocess.Popen(['hmmpfam2', '-E', str(max_eval), hmm_dir, '-'], stdin=stdin, stdout=subprocess.PIPE) as process:
        output, error = process.communicate(input=input_bytes)
        return StringIO(output.decode())



    # with open(hmm_out, 'w') as out:
    #     command = ['hmmpfam2', '-E', str(max_eval), hmm_dir, '-']
    #     subprocess.call(command, stdin=fasta_sequence, stdout=out)
