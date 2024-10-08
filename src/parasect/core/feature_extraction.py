import os
import subprocess
import warnings
from typing import Union, List

from Bio import SearchIO
from Bio.SearchIO._model import HSP
from pikachu.fingerprinting.ecfp_4 import ECFP
from pikachu.general import read_smiles, Structure

from parasect.core.parsers import parse_morgan_fingerprint, parse_amino_acid_properties, read_fasta


ONE_HOT_CATEGORIES = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
SEPARATOR_1 = '|'
SEPARATOR_2 = '_'
SEPARATOR_3 = '-'
START_POSITION = 66
MAPPING_SUFFIX: str = 'mapping.txt'
FASTA_SUFFIX: str = 'renamed_fasta.txt'
REF_SEQUENCE = "BAA00406.1.A1"
VALID_CHARACTERS = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-"}



import parasect.data

ALIGNMENT_FILE = os.path.join(os.path.dirname(parasect.data.__file__), 'structure_alignment.fasta')
PROPERTIES_FILE: str = os.path.join(os.path.dirname(parasect.data.__file__), 'physicochemical_properties.txt')
PROPERTIES: dict[str, list[float]] = parse_amino_acid_properties(PROPERTIES_FILE)
HMM2_FILE = os.path.join(os.path.dirname(parasect.data.__file__), 'AMP-binding_hmmer2.hmm')
REF_SEQ_FILE = os.path.join(os.path.dirname(parasect.data.__file__), 'reference_sequence.fasta')


APOSITION_FILE = os.path.join(os.path.dirname(parasect.data.__file__), 'stachelhaus.txt')
APOSITION_FILE_34 = os.path.join(os.path.dirname(parasect.data.__file__), 'active_site.txt')
APOSITION_FILE_HMM2 = os.path.join(os.path.dirname(parasect.data.__file__), 'stachelhaus_hmm2.txt')
APOSITION_FILE_34_HMM2 = os.path.join(os.path.dirname(parasect.data.__file__), 'active_site_hmm2.txt')



def bitvector_from_smiles(smiles_string: str, bitvector_file: str) -> list[int]:
    """
    Return a bitvector from a SMILES string

    Parameters
    ----------
    smiles_string: str, must be a valid SMILES string
    bitvector_file: path to file containing bitvectors on which paras and parasect have been trained

    Returns
    -------
    vector: list of int, with each int either 0 or 1, 0 denoting absence and 1 denoting prescence of bitvector. Vector
        is given in the same order as in bitvector_file

    """
    structure: Structure = read_smiles(smiles_string)

    ecfp: ECFP = ECFP(structure)
    fingerprint: set[int] = ecfp.fingerprint

    with open(bitvector_file, 'r') as bitvectors:
        substructure_hashes = list(map(int, bitvectors.readline().split('\t')[2:]))

    vector: list[int] = []

    for substructure_hash in substructure_hashes:
        if substructure_hash in fingerprint:
            vector.append(1)
        else:
            vector.append(0)

    return vector


def bitvectors_from_substrate_names(substrate_names, fingerprint_file):
    """
    Return substrate names and fingerprints for all substrate names for which a fingerprint could be found

    Parameters
    ----------
    substrate_names: list of [str, ->], with each str a substrate name. If the substrate name is not in the fingerprint
        file, a warning will be printed.
    fingerprint_file: str, path to file containing precomputed fingerprints, with one substrate per row and one
        substructure per column

    Returns
    -------
    substrates: list of [str, ->], substrate names which were present in the fingerprint file
    fingerprints: list of [[int, ->], ->], with each list of integers a fingerprint. The index of each fingerprint
        matches the index of the corresponding substrate name in substrates.

    """
    substrate_to_fingerprint: dict[str, list[int]] = parse_morgan_fingerprint(fingerprint_file)
    substrates: list[str] = []
    fingerprints: list[list[int]] = []
    for substrate in substrate_names:
        if substrate in substrate_to_fingerprint:
            substrates.append(substrate)
            fingerprints.append(substrate_to_fingerprint[substrate])
        else:
            warnings.warn(f"Could not find a fingerprint for substrate name {substrate}. Excluded from analysis.")
    return substrates, fingerprints


def run_hmmpfam2(hmm_dir, fasta_file, hmm_out):
    """
    Run hmmscan from command line

    Input:
    hmm_dir: str, dir of .hmm file containing the HMMs to be used in the scan
    fasta_dir: str, dir of .fasta file containing the sequences to be scanned
    out_dir: str, file location containing results of hmmscan

    """

    with open(hmm_out, 'w') as out:
        command = ['hmm2pfam', hmm_dir, fasta_file]
        # command = ['hmmpfam2', hmm_dir, fasta_file]
        # command = [hmmer_path, hmm_dir, fasta_file]
        subprocess.call(command, stdout=out)



def run_muscle(in_file, alignment_file, out_file):

    command = ['muscle', '-quiet', '-profile',  '-in1', alignment_file, '-in2', in_file, '-out', out_file]

    subprocess.check_call(command)


def align_adomain(domain_name, domain_sequence, alignment_file, temp_dir):
    # Run muscle and collect sequence positions from file

    temp_in = os.path.join(temp_dir, 'temp_in_alignment.fasta')
    temp_out = os.path.join(temp_dir, 'temp_out_alignment.fasta')

    with open(temp_in, 'w') as temp:
        temp.write(f'>{domain_name}\n{domain_sequence}')

    # Align sequence to all a domains in database with muscle
    run_muscle(temp_in, alignment_file, temp_out)
    id_to_alignment = read_fasta(temp_out)

    # Aligned sequence of domain
    aligned_domain = id_to_alignment[domain_name]

    # Aligned sequence of 1AMU reference sequence
    aligned_reference = id_to_alignment[REF_SEQUENCE]

    return aligned_domain, aligned_reference


class AdenylationDomain:
    def __init__(self, protein_name, domain_start, domain_end):
        self.protein_name = protein_name
        self.start = domain_start
        self.end = domain_end

        self.sequence = None
        self.signature = None
        self.extended_signature = None
        self.domain_nr = 0
        self.domain_id = None

    def set_domain_id(self, separator_1, separator_2, separator_3):
        assert self.protein_name
        assert self.domain_nr

        self.domain_id = f"{self.protein_name}{separator_1}domain{separator_2}{self.domain_nr}{separator_1}{self.start}{separator_3}{self.end}"

    def set_domain_signatures_hmm(self, hit_n_terminal, hit_c_terminal=None):
        """Extract (extended) signatures from adenylation domains using HMM"""

        signature_positions = HMM2_POSITIONS_SIGNATURE
        extended_signature_positions = HMM2_POSITIONS_EXTENDED_SIGNATURE
        position_k = HMM2_POSITION_K

        profile = hit_n_terminal.aln[1].seq
        query = hit_n_terminal.aln[0].seq
        offset = hit_n_terminal.hit_start

        signature = get_reference_positions_hmm(query, profile, [p - offset for p in signature_positions])
        if signature and all([char in VALID_CHARACTERS for char in signature]):
            self.signature = signature

        lysine = None

        if hit_c_terminal:
            profile_c = hit_c_terminal.aln[1].seq
            query_c = hit_c_terminal.aln[0].seq
            offset_c = hit_c_terminal.hit_start
            lysine = get_reference_positions_hmm(query_c, profile_c, [p - offset_c for p in position_k])

        if self.signature:
            if lysine and lysine in VALID_CHARACTERS and lysine != '-':
                self.signature += lysine
            else:
                self.signature += "K"

        extended_signature = get_reference_positions_hmm(query, profile, [p - offset for p in extended_signature_positions])
        if extended_signature and all([char in VALID_CHARACTERS for char in extended_signature]):
            self.extended_signature = extended_signature

    def set_domain_signatures_profile(self, temp_dir):
        """Extract (extended) signatures from adenylation domains using profile alignment"""

        if not self.sequence:
            raise Exception("Sequence needs to be defined first")

        seq_id = "DOMAIN_TO_QUERY"

        aligned_domain, aligned_reference = align_adomain(seq_id, self.sequence, ALIGNMENT_FILE, temp_dir)

        aligned_positions_signature = get_reference_positions(POSITIONS_SIGNATURE, aligned_reference)
        aligned_positions_extended_signature = get_reference_positions(POSITIONS_EXTENDED_SIGNATURE, aligned_reference)

        signature = []
        for position in aligned_positions_signature:
            signature.append(aligned_domain[position])
        self.signature = ''.join(signature)

        extended_signature = []
        for position in aligned_positions_extended_signature:
            extended_signature.append(aligned_domain[position])
        self.extended_signature = ''.join(extended_signature)

    def set_domain_number(self, domain_nr):
        self.domain_nr = domain_nr

    def set_sequence(self, sequence):
        self.sequence = sequence

    def write_sequence(self, fasta_file, sequence_type='full'):
        assert self.domain_id
        if sequence_type == 'full':
            fasta_file.write(f'>{self.domain_id}\n{self.sequence}\n')
        elif sequence_type == 'extended_signature':
            fasta_file.write(f'>{self.domain_id}\n{self.extended_signature}\n')
        elif sequence_type == 'signature':
            fasta_file.write(f'>{self.domain_id}\n{self.signature}\n')
        else:
            raise ValueError(f"Only accepted sequence types are 'full', 'signature', and 'extended_signature'. Got {sequence_type}")


def rename_sequences(fasta_file: str, out_dir: str) -> tuple[str, str]:
    """

    Rename sequences before running hmmscan, and return file paths of mapping and new fasta file

    Parameters
    ----------
    fasta_file: str, path to input fasta file
    out_dir: str, path to output directory

    Returns
    -------
    mapping_file: str, path to mapping file which maps renamed sequence IDs to the original sequence IDs
    new_fasta_file: str, path to output fasta file

    """
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    mapping_file: str = os.path.join(out_dir, MAPPING_SUFFIX)
    new_fasta_file: str = os.path.join(out_dir, FASTA_SUFFIX)

    id_to_seq: dict[str, str] = read_fasta(fasta_file)
    counter: int = 0
    with open(new_fasta_file, 'w') as new_fasta:
        with open(mapping_file, 'w') as mapping:
            for seq_id, seq in id_to_seq.items():
                counter += 1
                mapping.write(f"{counter}\t{seq_id}\n")
                new_fasta.write(f">{counter}\n{seq}\n")

    return mapping_file, new_fasta_file


def parse_mapping(mapping_file):
    """
    Return a dictionary of renamed sequence IDs to original sequence IDs

    Parameters
    ----------
    mapping_file: str, path to mapping file which maps renamed sequence IDs to the original sequence IDs

    Returns
    -------
    new_to_original: dict of {new_id: old_id, ->}, with new_id and old_id strings

    """
    new_to_original = {}
    with open(mapping_file, 'r') as mapping:
        for line in mapping:
            line = line.strip()
            line_segments = line.split('\t')
            new = line_segments[0]
            original = '\t'.join(line_segments[1:])
            new_to_original[new] = original

    return new_to_original


def reverse_renaming(adenylation_domains: list["AdenylationDomain"], mapping_file: str) -> None:
    """
    Reverses the renaming of sequences within adenylation domain instances

    Parameters
    ----------

    adenylation_domains: list of [domain, ->], with each domain an AdenylationDomain instance
    mapping_file: str, path to mapping file which maps renamed sequence IDs to the original sequence IDs

    """
    new_to_original: dict[str, str] = parse_mapping(mapping_file)
    for domain in adenylation_domains:
        domain.protein_name = new_to_original[domain.protein_name]


def domains_to_features(domains: list["AdenylationDomain"], signature: str = "extended", one_hot: bool = False):
    """
    Return a list of sequence ids and a list of corresponding feature vectors from a list of AdenylationDomain instances

    Parameters
    ----------
    domains: list of [domain, ->], with each domain an AdenylationDomain instance
    signature: str, type of signature used for featurisation. Must be 'extended' or 'short'
    one_hot: bool, use one-hot encoding if True, the 15 NRPSPredictor features otherwise

    Returns
    -------
    sequence_ids: list of [sequence_id, ->], with each sequence_id str
    feature_vectors: list of [[feature, ->], ->], with each feature float or int. Index of feature list corresponds to
        index of corresponding sequence_id in sequence_ids

    """
    sequence_ids: list[str] = []
    feature_vectors: list[Union[list[float], list[int]]] = []

    for domain in domains:
        sequence_ids.append(domain.domain_id)

        if signature == "extended":
            if not one_hot:
                feature_vector: list[float] = get_sequence_features(domain.extended_signature)
            else:
                feature_vector: list[int] = one_hot_encoding(domain.extended_signature)

        elif signature == "short":
            if not one_hot:
                feature_vector: list[float] = get_sequence_features(domain.signature)
            else:
                feature_vector: list[int] = one_hot_encoding(domain.signature)
        else:
            raise ValueError(f"Expected 'extended' or 'short'. Got {signature}.")

        feature_vectors.append(feature_vector)

    return sequence_ids, feature_vectors


def get_sequence_features(sequence: str) -> list[float]:
    """
    Return a feature list of NRPSPredictor features from a sequence

    Parameters
    ----------
    sequence: str, amino acid sequence

    Returns
    -------
    features: list of [feature, ->], with each feature float

    """
    features: list[float] = []

    for aa in sequence:
        properties = PROPERTIES[aa]
        features.extend(properties)

    return features


def one_hot_encoding(sequence: str) -> list[int]:
    """
    Return one-hot encoding of a sequence

    Parameters
    ----------
    sequence: str, amino acid sequence

    Returns
    -------
    feature_vector: list of [feature, ->], with each feature int

    """
    feature_vector: list[int] = []

    for res in sequence:
        assert (res in ONE_HOT_CATEGORIES or res == '-')
        for amino_acid in ONE_HOT_CATEGORIES:

            if res == amino_acid:
                feature_vector.append(1)
            else:
                feature_vector.append(0)

    return feature_vector


def parse_hmm2_results(hmm_results):
    """
    Return dictionary of domain identifier to Biopython HSP instance

    Parameters
    ----------
    hmm_results: str, path to hmmpfam2 output file (hmmer-2)

    Returns
    -------
    id_to_hit: dict of {domain_id: HSP, ->}, with domain_id str and HSP a Biopython HSP instance containing a HMM hit
        between an Hmm2 adenylation domain HMM and a query sequence

    """
    id_to_hit: dict[str, HSP] = {}

    for result in SearchIO.parse(hmm_results, 'hmmer2-text'):
        for hsp in result.hsps:
            if hsp.bitscore > 20:
                if hsp.hit_id == 'AMP-binding' or hsp.hit_id == 'AMP-binding_C':
                    header = make_domain_id(result.id, hsp.hit_id, hsp.query_start, hsp.query_end)
                    id_to_hit[header] = hsp

    return id_to_hit


def hits_to_domains(id_to_hit, fasta_file, temp_dir, profile=False, verbose=False):
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
                a_domain.set_domain_signatures_profile(temp_dir)
                if a_domain.sequence and a_domain.extended_signature and a_domain.signature and a_domain.domain_nr:
                    filtered_a_domains.append(a_domain)

        filtered_a_domains.sort(key=lambda x: (x.protein_name, x.start))

    return filtered_a_domains


def domains_from_fasta(fasta_in, temp_dir, job_name='paras_run', profile=False, verbose=False):
    """
    Extract adomains from fasta file and write them to out_dir

    Input:
    fasta_in: str, file location of .fasta file containing aa sequences
    job_name: str, name of job

    Output:
    fasta_out: str, file location of .fasta file containing detected adomains

    """
    hmm_out = os.path.join(temp_dir, f'{job_name}.hmm_result')

    if verbose:
        print("Running HMM..")
    run_hmmpfam2(HMM2_FILE, fasta_in, hmm_out)
    if verbose:
        print("Parsing hits..")
    id_to_hit = parse_hmm2_results(hmm_out)

    if profile:
        if verbose:
            print("Processing hits (profile alignment-based active site extraction)..")
        a_domains = hits_to_domains(id_to_hit, fasta_in, temp_dir, profile=True, verbose=verbose)
    else:
        if verbose:
            print("Processing hits (hmm-based active site extraction)..")
        a_domains = hits_to_domains(id_to_hit, fasta_in, temp_dir, verbose=verbose)

    return a_domains


if __name__ == "__main__":
    domains_from_fasta(REF_SEQ_FILE)


def read_positions(position_file: str, start_position: int) -> List[int]:
    """
    Return positions from a tab-separated file. Positions are relative to start_position.

    Input:
    position_file: str, the path to tab-separated file containing the positions
    start_position: int, a relative start position to adjust all positions by

    Output:
    positions: list of [int, ->], one for each position found in the file
    """

    with open(position_file, 'r') as position_data:
        text = position_data.read().strip()
        positions = []
        for i in text.split("\t"):
            positions.append(int(i) - start_position)
    return positions


def get_reference_positions(positions: List[int], aligned_reference: str) -> List[int]:
    """
    Adjusts a list of positions to account for gaps in the reference sequence

    Input:
    positions: list of [int, ->], with integers representing positions of interest in
        the reference sequence
    aligned_reference: the (aligned) reference sequence

    Output:
    pos_list: list of [int, ->], a new list of positions, each >= the original position
    """
    pos_list = []
    position = 0
    for i, aa in enumerate(aligned_reference):
        if aa != "-":
            if position in positions:
                pos_list.append(i)
            position += 1
    return pos_list


# From antiSMASH 7.0
def get_reference_positions_hmm(query, reference, ref_positions):
    """ Extracts the given positions from a query alignment. The positions are
        adjusted to account for any gaps in the reference sequence.

        Arguments:
            query: the aligned query
            reference: the aligned reference
            ref_positions: the positions of interest in the unaligned reference

        Returns:
            a string containing the sequence elements at the adjusted reference
            positions or None if the reference is too short for some reason
    """
    # adjust position of interest to account for gaps in the ref sequence alignment
    positions = []
    position_skipping_gaps = 0
    for i, amino in enumerate(reference):
        if amino in "-.":
            continue
        if position_skipping_gaps in ref_positions:
            positions.append(i)
        position_skipping_gaps += 1
    if len(positions) != len(ref_positions):
        return None
    # extract positions from query sequence
    return "".join([query[i] for i in positions])


POSITIONS_SIGNATURE = read_positions(APOSITION_FILE, START_POSITION)
POSITIONS_EXTENDED_SIGNATURE = read_positions(APOSITION_FILE_34, START_POSITION)

HMM2_POSITIONS_SIGNATURE = read_positions(APOSITION_FILE_HMM2, 0)
HMM2_POSITIONS_EXTENDED_SIGNATURE = read_positions(APOSITION_FILE_34_HMM2, 0)


# Hmmer2:
HMM2_POSITION_K = [36]


def parse_domain_id(fasta_id):
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
    seq_id, hit_id, hit_location = fasta_id.split(SEPARATOR_1)
    hit_start, hit_end = hit_location.split('-')
    hit_start = int(hit_start)
    hit_end = int(hit_end)
    return seq_id, hit_id, hit_start, hit_end


def make_domain_id(seq_id, hit_id, start, end):
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
    return f'{seq_id}{SEPARATOR_1}{hit_id}{SEPARATOR_1}{start}-{end}'
