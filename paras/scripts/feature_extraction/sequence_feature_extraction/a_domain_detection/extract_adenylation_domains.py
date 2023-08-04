#!/usr/bin/env python

import subprocess
import os

from Bio import SearchIO

from paras.scripts.parsers.clear_temp import clear_temp
from paras.scripts.parsers.fasta import read_fasta

import paras.data.temp
import paras.data.sequence_data.hmm


HMM_FILE = os.path.join(os.path.dirname(paras.data.sequence_data.hmm.__file__), 'AMP-binding_full.hmm')
TEMP_DIR = os.path.dirname(paras.data.temp.__file__)


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


def remove_insertions(seq):
    """
    Remove insertion states from a sequence matched to an HMM

    Input:
    seq: str, amino acid sequence where uppercase characters are match states
        and lowercase characters are insertion states

    Output:
    new_seq: str, amino acid sequence with all insertion states (lowercase
        characters) removed.
    """
    new_seq = []
    for character in seq:
        if not character.islower():
            new_seq.append(character)

    new_seq = ''.join(new_seq)
    return new_seq


def make_header(ID, hit_id, start, end):
    """
    Return header (str) for .fasta output file

    Input:
    ID: str, sequence identifier
    hit_id: str, HMM the sequence matched to
    start: int, first aa in the query sequence that the HMM matches to
    end: int, last aa in the query sequence that the HMM matches to

    Output:
    str, sequence properties separated by '|'
    """
    return f'{ID}|===|{hit_id}|===|{start}-{end}'


def remove_gaps(sequence):
    """
    Remove gaps from sequence

    Input:
    sequence: str, amino acid sequence

    Output:
    str, same sequence as input with all '-' characters removed
    """
    new_sequence = []
    for character in sequence:
        if character != '-':
            new_sequence.append(character.upper())
    return ''.join(new_sequence)


def parse_hmm_results(hmm_results, fasta_out):
    """
    Write sequences that match to AMP-binding domains to .fasta file

    Input:
    hmm_results: str, file location containing results of hmmscan
    fasta_out: str, file location for writing AMP-binding domains
    """
    fasta_file = open(fasta_out, 'w')
    for result in SearchIO.parse(hmm_results, 'hmmer3-text'):
        for hsp in result.hsps:
            if hsp.evalue < 0.00001:
                if hsp.hit_id == 'AMP-binding' or hsp.hit_id == 'AMP-binding_C':
                    seq = remove_gaps(hsp.query.seq)
                    header = make_header(result.id, hsp.hit_id, hsp.query_start, hsp.query_end)
                    fasta_file.write(">%s\n%s\n" % (header, seq))

    fasta_file.close()


def parse_fasta_id(fasta_id):
    """
    Return id, domain type and domain location from common id

    Input:
    fasta_id: str

    Output:
    id: str, sequence id
    hit_id: str, domain id
    hit_start: int, start position of domain in protein
    hit_end: int, end position of domain in protei
    """
    id, hit_id, hit_location = fasta_id.split('|===|')
    hit_start, hit_end = hit_location.split('-')
    hit_start = int(hit_start)
    hit_end = int(hit_end)
    return id, hit_id, hit_start, hit_end


def merge_adomains(original_fasta_dir, fasta_dir, new_fasta_dir):
    """
    Find whole AMP-binding domains from hits to N- and C-terminal domains

    Consecutive hits to an AMP-binding domain and then an AMP-binding_C domain
    within 30 amino acids of each other are joined into a single AMP-binding
    domain. Standalone AMP-binding_C domains will be ignored, while standalone
    AMP-binding domains will still be reported.

    Input:
    fasta_dir: str, file location of .fasta file containing AMP-binding and
        AMP-binding_C domains
    original_fasta_dir: str, file location of .fasta file from which
        AMP-binding domains were extracted
    new_fasta_dir: str, file location of .fasta file to which whole
        AMP-binding domains will be written

    """
    fasta = read_fasta(fasta_dir)
    hits_by_seq_id = {}
    for id in fasta:
        seq_id, hit_id, hit_start, hit_end = parse_fasta_id(id)
        if not seq_id in hits_by_seq_id:
            hits_by_seq_id[seq_id] = []

        hits_by_seq_id[seq_id].append((hit_id, hit_start, hit_end))

    amp_nc_pairs = []

    for seq_id, hits in hits_by_seq_id.items():
        for hit_1 in hits:
            if hit_1[0] == 'AMP-binding':
                match_found = False
                for hit_2 in hits:
                    if hit_2[0] == 'AMP-binding_C':
                        if (hit_2[1] > hit_1[2]) and (hit_2[1] - hit_1[2] < 30):
                            match_found = True
                            amp_nc_pairs.append((seq_id, hit_1[1], hit_2[2]))

                if not match_found:
                    amp_nc_pairs.append((seq_id, hit_1[1], hit_1[2]))



    original_fasta = read_fasta(original_fasta_dir)
    amp_nc_pairs = sorted(amp_nc_pairs, key=lambda x: x[1])

    with open(new_fasta_dir, 'w') as new_fasta:
        for seq_id, sequence in original_fasta.items():
            counter = 0
            for amp_nc_pair in amp_nc_pairs:
                pair_id, start, end = amp_nc_pair
                if seq_id == pair_id:

                    amp_sequence = sequence[start:end]
                    if len(amp_sequence) > 280:
                        counter += 1
                        header = f'{seq_id}__{counter}__{start + 1}-{end}'

                        new_fasta.write(f'>{header}\n{amp_sequence}\n')


def find_adomains(fasta_in, out_dir):
    """
    Extract adomains from fasta file and write them to out_dir

    Input:
    fasta_in: str, file location of .fasta file containing aa sequences
    out_dir: str, directory where output file will be stored

    Output:
    fasta_out: str, file location of .fasta file containing detected adomains

    """
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    # Remove extension from file string
    file_label = '.'.join(fasta_in.split('.')[:-1])
    file_label = file_label.split(r'/')[-1]

    hmm_out = os.path.join(TEMP_DIR, f'{file_label}_.hmm_result')
    fasta_subdomains_out = os.path.join(TEMP_DIR, f'{file_label}_subdomains.fasta')
    fasta_out = os.path.join(out_dir, f'{file_label}_adomains.fasta')

    run_hmmscan(HMM_FILE, fasta_in, hmm_out)
    parse_hmm_results(hmm_out, fasta_subdomains_out)
    merge_adomains(fasta_in, fasta_subdomains_out, fasta_out)
    clear_temp()

    return fasta_out

