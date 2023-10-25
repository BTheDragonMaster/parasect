#!/usr/bin/env python

import os

from Bio import SearchIO

from paras.scripts.parsers.fasta import read_fasta
from paras.scripts.data_processing.temp import TEMP_DIR
from paras.scripts.feature_extraction.sequence_feature_extraction.hmm.run_hmmscan import HMM2_FILE, \
    run_hmmpfam2
from paras.scripts.feature_extraction.sequence_feature_extraction.sequence_labels import parse_domain_id, \
    make_domain_id
from paras.scripts.feature_extraction.sequence_feature_extraction.adenylation_domain import AdenylationDomain

import paras.data.sequence_data.sequences

REF_SEQ_FILE = os.path.join(os.path.dirname(paras.data.sequence_data.sequences.__file__), 'reference_sequence.fasta')


def parse_hmm3_results(hmm_results):
    """
    Write sequences that match to AMP-binding domains to .fasta file

    Input:
    hmm_results: str, file location containing results of hmmscan
    fasta_out: str, file location for writing AMP-binding domains
    """
    id_to_hit = {}

    for result in SearchIO.parse(hmm_results, 'hmmer3-text'):
        for hsp in result.hsps:
            if hsp.bitscore > 20:
                if hsp.hit_id == 'AMP-binding' or hsp.hit_id == 'AMP-binding_C':
                    header = make_domain_id(result.id, hsp.hit_id, hsp.query_start, hsp.query_end)
                    id_to_hit[header] = hsp

    return id_to_hit


def parse_hmm2_results(hmm_results):
    """
    Write sequences that match to AMP-binding domains to .fasta file

    Input:
    hmm_results: str, file location containing results of hmmscan
    fasta_out: str, file location for writing AMP-binding domains
    """
    id_to_hit = {}

    for result in SearchIO.parse(hmm_results, 'hmmer2-text'):
        for hsp in result.hsps:
            if hsp.bitscore > 20:
                if hsp.hit_id == 'AMP-binding' or hsp.hit_id == 'AMP-binding_C':
                    header = make_domain_id(result.id, hsp.hit_id, hsp.query_start, hsp.query_end)
                    id_to_hit[header] = hsp

    return id_to_hit


def hits_to_domains(id_to_hit, fasta_file, profile=False, verbose=False):
    hits_by_seq_id = {}
    if verbose:
        print("\tSorting hits by sequence..")

    for hit_key, hit in id_to_hit.items():
        seq_id, hit_id, hit_start, hit_end = parse_domain_id(hit_key)
        if seq_id not in hits_by_seq_id:
            hits_by_seq_id[seq_id] = []

        hits_by_seq_id[seq_id].append((hit_id, hit_start, hit_end, hit_key))

    if verbose:
        print("\tExtracting domain signatures..")

    counter = 0

    seq_id_to_domains = {}

    for seq_id, hits in hits_by_seq_id.items():
        counter += 1
        for hit_id_1, hit_start_1, hit_end_1, hit_key_1 in hits:
            if hit_id_1 == 'AMP-binding':
                if seq_id not in seq_id_to_domains:
                    seq_id_to_domains[seq_id] = []
                match_found = False
                for hit_id_2, hit_start_2, hit_end_2, hit_key_2 in hits:
                    if hit_id_2 == 'AMP-binding_C':
                        # Todo: check that 200 is a good cutoff score

                        if (hit_start_2 > hit_end_1) and (hit_start_2 - hit_end_1 < 200):
                            a_domain = AdenylationDomain(seq_id, hit_start_1, hit_end_2)
                            if not profile:
                                a_domain.set_domain_signatures_hmm(id_to_hit[hit_key_1], id_to_hit[hit_key_2])
                            seq_id_to_domains[seq_id].append(a_domain)
                            match_found = True
                            break

                if not match_found:
                    a_domain = AdenylationDomain(seq_id, hit_start_1, hit_end_1)
                    a_domain.set_domain_signatures_hmm(id_to_hit[hit_key_1])
                    seq_id_to_domains[seq_id].append(a_domain)

        if verbose and counter % 1000 == 0:
            print(f"\t\tProcessed {counter} proteins.")

    for seq_id, domains in seq_id_to_domains.items():
        domains.sort(key=lambda x: x.start)

    fasta = read_fasta(fasta_file)

    if verbose:
        print("\tSetting domain numbers..")

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

    if not profile:
        if verbose:
            print("\tFiltering domains..")

        filtered_a_domains = []

        for seq_id, a_domains in seq_id_to_domains.items():
            for a_domain in a_domains:
                if a_domain.sequence and a_domain.extended_signature and a_domain.signature and a_domain.domain_nr:
                    filtered_a_domains.append(a_domain)

        if verbose:
            print("\tSorting domains by protein name..")

        filtered_a_domains.sort(key=lambda x: (x.protein_name, x.start))

    else:
        if verbose:
            print("\tExtracting domain signatures with profile alignment..")

        filtered_a_domains = []
        for seq_id, a_domains in seq_id_to_domains.items():
            for a_domain in a_domains:
                a_domain.set_domain_signatures_profile()
                if a_domain.sequence and a_domain.extended_signature and a_domain.signature and a_domain.domain_nr:
                    filtered_a_domains.append(a_domain)

        filtered_a_domains.sort(key=lambda x: (x.protein_name, x.start))

    return filtered_a_domains


def domains_from_fasta(fasta_in, job_name='paras_run', profile=False, verbose=False):
    """
    Extract adomains from fasta file and write them to out_dir

    Input:
    fasta_in: str, file location of .fasta file containing aa sequences
    job_name: str, name of job

    Output:
    fasta_out: str, file location of .fasta file containing detected adomains

    """

    hmm_out = os.path.join(TEMP_DIR, f'{job_name}.hmm_result')
    if verbose:
        print("Running HMM..")
    run_hmmpfam2(HMM2_FILE, fasta_in, hmm_out)
    if verbose:
        print("Parsing hits..")
    id_to_hit = parse_hmm2_results(hmm_out)

    if profile:
        if verbose:
            print("Processing hits (profile alignment-based active site extraction)..")
        a_domains = hits_to_domains(id_to_hit, fasta_in, profile=True, verbose=verbose)
    else:
        if verbose:
            print("Processing hits (hmm-based active site extraction)..")
        a_domains = hits_to_domains(id_to_hit, fasta_in, verbose=verbose)

    return a_domains


if __name__ == "__main__":
    domains_from_fasta(REF_SEQ_FILE)
