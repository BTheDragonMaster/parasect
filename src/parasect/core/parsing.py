# -*- coding: utf-8 -*-

"""Parsers module for PARASECT."""

import os
import itertools
from collections import OrderedDict
from typing import Dict, List, Union

from Bio import SeqIO


class Tabular:
    """Class for reading and storing tabular data files."""

    def __init__(self, path_in: str, separator: str = "\t") -> None:
        """Initializes the Tabular class.
        
        :param path_in: Path to tabular file.
        :type path_in: str
        :param separator: Separator used in the tabular file (default is tab).
        :type separator: str
        :raises FileNotFoundError: If the file at the specified path does not exist.
        :raises ValueError: If the number of columns in a row does not match the 
            length of the number of columns in the header.
        :raises ValueError: If there is a duplicate row ID in the data.
        """
        self.index_of_column_with_id = 0  # default index of column with ID
        self.column_names = []
        self.rows = OrderedDict()

        # check if the file exists
        if not os.path.exists(path_in):
            raise FileNotFoundError(f"file not found: {path_in}")

        # read the tabular file
        with open(path_in, "r") as fo:
            self.column_names = [n.strip() for n in fo.readline().strip().split(separator)]

            for line_idx, line in enumerate(fo):
                # split the line into a list of values
                row = [v.strip() for v in line.strip().split(separator)]

                # check if the row has the same number of columns as the header
                if len(row) != len(self.column_names):
                    msg = f"row {line_idx + 2} has a different number of columns than the header"  # noqa: E501
                    raise ValueError(msg)

                # parse out the row ID
                row_id = row[self.index_of_column_with_id]

                # check for duplicate row IDs
                if row_id in self.rows:
                    msg = f"duplicate row ID when reading {path_in}: {row_id}"
                    raise ValueError(msg)

                # if the row ID is unique, add the row to the data dictionary
                self.rows[row_id] = OrderedDict()

                for value_idx, value in enumerate(row):
                    column_name = self.column_names[value_idx]
                    self.rows[row_id][column_name] = value

    def get_column_values(self, column_name: str) -> List[Union[int, float, str]]:
        """Returns a list of values from a specified column.
        
        :param column_name: Name of the column.
        :type column_name: str
        :return: List of values from the specified column.
        :rtype: List[Union[int, float, str]]
        :raises KeyError: If the specified column name is not found in the data.
        """
        # check if the column name is in the column names
        if column_name not in self.column_names:
            raise KeyError(f"cannot find category {column_name} in data")

        # get the values from the specified column
        column_values = []
        for row_id in self.rows:
            column_values.append(self.get_row_value(row_id, column_name))

        return column_values

    def get_row_values(self, row_id: str) -> List[Union[int, float, str]]:
        """Returns a row of values from the data.
        
        :param row_id: ID of the row.
        :type row_id: str
        :return: Row of values from the data.
        :rtype: List[Union[int, float, str]]
        :raises KeyError: If the specified row ID is not found in the data.
        """
        # check if the row ID is in the data
        if row_id not in self.rows:
            raise KeyError(f"cannot find data ID {row_id} in data")

        # get the values from the specified row
        row_values = []
        for category in self.column_names:
            row_values.append(self.get_row_value(row_id, category))

        return row_values

    def get_row_value(self, row_id: str, column_name: str) -> Union[int, float, str]:
        """Returns a value from a specified row and column.
        
        :param row_id: ID of the row.
        :type row_id: str
        :param column_name: Name of the column.
        :type column_name: str
        :return: Value from the specified row and column.
        :rtype: Union[int, float, str]
        :raises KeyError: If the specified row ID or column name is not found in the data.
        """
        return self.rows[row_id][column_name]

    def write_table(self, path_out: str, separator: str = "\t") -> None:
        """Writes the tabular data to a file.

        :param path_out: Path to output tabular file.
        :type path_out: str
        :param separator: Separator to use in the output file (default is tab).
        :type separator: str
        """
        with open(path_out, "w") as fo:

            # write the header
            fo.write(separator.join(self.column_names) + "\n")

            # write the data
            for row_id in self.rows:
                row_values = self.rows[row_id]
                fo.write(separator.join(row_values.values()) + "\n")


def parse_substrate_list(path_in: str) -> List[str]:
    """Parses a list of substrate names from a file.
    
    :param path_in: Path to input file.
    :type path_in: str
    :return: List of substrates.
    :rtype: List[str]
    :raises FileNotFoundError: If the file at the specified path does not exist.
    """
    # check if the file exists
    if not os.path.exists(path_in):
        raise FileNotFoundError(f"substrate list file not found: {path_in}")

    # initialize the list of substrates
    substrate_names = []

    # read the file
    with open(path_in, "r") as fo:
        for line in fo:
            substrate_name = line.strip()

            # check if the line is not empty
            if substrate_name:
                substrate_names.append(substrate_name)

    return substrate_names


def parse_amino_acid_properties_file(path_in: str) -> Dict[str, List[float]]:
    """Parses a dictionary of amino acid to property vector from a file.
    
    :param path_in: Path to input file.
    :type path_in: str
    :return: Dictionary of amino acid to property vector.
    :rtype: Dict[str, List[float]]
    :raises FileNotFoundError: If the file at the specified path does not exist.
    """
    # check if the file exists
    if not os.path.exists(path_in):
        raise FileNotFoundError(f"amino acid properties file not found: {path_in}")

    # read the tabular file
    data = Tabular(path_in)

    # parse the amino acid properties from the file
    amino_acid_properties = {}
    for row_id in data.rows:

        # get the amino acid ID and property vector
        amino_acid_id = data.get_row_value(row_id, column_name="AA")
        amino_acid_props = [float(v) for v in data.get_row_values(row_id)[1:]]
        
        # check if the property vector is of length 15
        if len(amino_acid_props) != 15:
            raise ValueError(f"amino acid property vector for {amino_acid_id} is not of length 15")  # noqa: E501

        # add the amino acid ID and property vector to the dictionary
        amino_acid_properties[amino_acid_id] = amino_acid_props

    return amino_acid_properties


def parse_morgan_fingerprint_file(path_in: str) -> Dict[str, List[int]]:
    """Parses a dictionary of molecule name to Morgan fingerprint from a file.
    
    :param path_in: Path to input file.
    :type path_in: str
    :return: Dictionary of molecule name to Morgan fingerprint.
    :rtype: Dict[str, List[int]]
    :raises FileNotFoundError: If the file at the specified path does not exist.
    """
    # check if the file exists
    if not os.path.exists(path_in):
        raise FileNotFoundError(f"morgan fingerprint file not found: {path_in}")

    # read the tabular file
    data = Tabular(path_in)

    # parse the Morgan fingerprints from the file
    amino_acid_fingerprints = OrderedDict()
    for row_id in data.rows:
        # skip the molecule name and SMILESMILES
        fingerprint = [int(v) for v in data.get_row_values(row_id)[2:]] 

        # the row_id is the molecule name
        amino_acid_fingerprints[row_id] = fingerprint
    
    # print(amino_acid_fingerprints)
    return amino_acid_fingerprints


def read_fasta_file(path_in: str) -> Dict[str, str]:
    """Reads fasta file and parses it into a dictionary format.

    :param path_in: Path to fasta file.
    :type path_in: str
    :return: Dictionary of fasta sequences, where the key is the sequence header
        and the value is the sequence.
    :rtype: Dict[str, str]
    :raises FileNotFoundError: If the file at the specified path does not exist.
    """
    # check if the file exists
    if not os.path.exists(path_in):
        raise FileNotFoundError(f"File not found: {path_in}")

    # initialize the dictionary
    fasta_dict = {}
    
    with open(path_in, "r") as fo:

        # initialize the header and sequence list
        header = None
        sequence = []

        for line in fo:
            line = line.strip()

            if line.startswith(">"):
                if sequence:
                    # if the sequence list is not empty, join the list into a string
                    # and add it to the dictionary with the header as the key
                    fasta_dict[header] = "".join(sequence)

                    header = line[1:]  # remove the ">" from the header
                    sequence = []  # reset the sequence list
                else:
                    header = line[1:]

            else:
                sequence.append(line)
        
        # if the header is not None, this means there is a sequence that has not 
        # yet been added to the dictionary
        if header is not None:
            fasta_dict[header] = "".join(sequence)

    return fasta_dict


def write_fasta_file(fasta_dict: Dict[str, str], path_out: str) -> None:
    """Writes a dictionary of fasta sequences to a file.

    :param fasta_dict: Dictionary of fasta sequences, where the key is the sequence 
        header and the value is the sequence.
    :type fasta_dict: Dict[str, str]
    :param path_out: Path to output fasta file.
    :type path_out: str
    """
    with open(path_out, "w") as fo:
    
        # iterate over the dictionary items
        for header, sequence in fasta_dict.items():
            fo.write(f">{header}\n{sequence}\n")


def parse_protein_sequences_from_genbank_file(path_in: str, path_out: str) -> None:
    """Parses protein sequences from a GenBank file and writes them to a fasta file.

    :param path_in: Path to input GenBank file.
    :type path_in: str
    :param path_out: Path to output fasta file.
    :type path_out: str
    :raises FileNotFoundError: If the file at the specified path does not exist.
    """
    # check if the file exists
    if not os.path.exists(path_in):
        raise FileNotFoundError(f"File not found: {path_in}")

    # initialize a counter for generating gene IDs
    counter = itertools.count()

    # initialize a dictionary to store the fasta sequences
    fasta_dict = {}

    # parse the GenBank file
    for record in SeqIO.parse(path_in, "genbank"):
        for feature in record.features:

            # check if the feature is a coding sequence (CDS)
            if feature.type == "CDS":

                # check if the feature has a translation
                if "translation" in feature.qualifiers:
                    sequence = feature.qualifiers["translation"]

                    # check if the feature has a protein ID, gene ID, or locus tag
                    # to use as the sequence ID. If not, generate a gene ID
                    if "protein_id" in feature.qualifiers:
                        seq_id = feature.qualifiers["protein_id"][0]
                    elif "gene_id" in feature.qualifiers:
                        seq_id = feature.qualifiers["gene_id"][0]
                    elif "locus_tag" in feature.qualifiers:
                        seq_id = feature.qualifiers["locus_tag"][0]
                    else:
                        seq_id = f"gene_{next(counter)}"

                    # add the sequence to the dictionary
                    fasta_dict[seq_id] = sequence

    # write the fasta sequences to a file
    write_fasta_file(fasta_dict=fasta_dict, path_out=path_out)



####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################

MAPPING_SUFFIX: str = "mapping.txt"
FASTA_SUFFIX: str = "renamed_fasta.txt"


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

    id_to_seq: dict[str, str] = read_fasta_file(fasta_file)
    counter: int = 0
    with open(new_fasta_file, "w") as new_fasta:
        with open(mapping_file, "w") as mapping:
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
    with open(mapping_file, "r") as mapping:
        for line in mapping:
            line = line.strip()
            line_segments = line.split("\t")
            new = line_segments[0]
            original = "\t".join(line_segments[1:])
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


def domains_from_fasta(fasta_in, temp_dir, job_name="paras_run", profile=False, verbose=False):
    """
    Extract adomains from fasta file and write them to out_dir

    Input:
    fasta_in: str, file location of .fasta file containing aa sequences
    job_name: str, name of job

    Output:
    fasta_out: str, file location of .fasta file containing detected adomains

    """
    hmm_out = os.path.join(temp_dir, f"{job_name}.hmm_result")

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


def write_results(results, out_file):
    with open(out_file, "w") as out:
        out.write(f"sequence_id")
        header_written = False
        for seq_id, probabilities_and_substrates in results.items():
            if not header_written:
                for i in range(len(probabilities_and_substrates)):
                    out.write(f"\tsubstrate_{i + 1}\tconfidence_score_{i + 1}")
                out.write("\n")
            header_written = True

            out.write(f"{seq_id}")
            for probability, substrate in probabilities_and_substrates:
                out.write(f"\t{substrate}\t{probability:.2f}")
            out.write("\n")


def get_top_n_aa_paras(amino_acid_classes, probabilities, n):
    probs_and_aa = []
    for i, probability in enumerate(probabilities):
        probs_and_aa.append((probability, amino_acid_classes[i]))

    probs_and_aa.sort(reverse=True)

    return probs_and_aa[:n]


def get_top_n_aa_parasect(seq_id, id_to_probabilities, n):
    probabilities = id_to_probabilities[seq_id]
    probabilities.sort(reverse=True)
    return probabilities[:n]


def get_domains(
    input_file,
    extraction_method,
    job_name,
    separator_1,
    separator_2,
    separator_3,
    verbose,
    file_type,
    temp_dir,
):
    assert extraction_method in [
        "hmm",
        "profile",
    ], f"Only supported extraction methods are 'hmm' or 'profile'. Got {extraction_method}."
    assert file_type in [
        "fasta",
        "gbk",
    ], f"Only supported file types are 'fasta' or 'gbk'. Got {file_type}."

    if file_type == "gbk":
        original_fasta = os.path.join(temp_dir, "proteins_from_genbank.fasta")
        parse_protein_sequences_from_genbank_file(input_file, original_fasta)  # TODO: ???
        mapping_file, renamed_fasta_file = rename_sequences(original_fasta, temp_dir)
    else:
        mapping_file, renamed_fasta_file = rename_sequences(input_file, temp_dir)

    if extraction_method == "profile":
        a_domains = domains_from_fasta(
            renamed_fasta_file, temp_dir=temp_dir, job_name=job_name, profile=True, verbose=verbose
        )
    elif extraction_method == "hmm":
        a_domains = domains_from_fasta(
            renamed_fasta_file, temp_dir=temp_dir, job_name=job_name, verbose=verbose
        )
    else:
        raise ValueError(
            f"Only supported extraction methods are 'hmm' or 'profile'. Got {extraction_method}."
        )

    reverse_renaming(a_domains, mapping_file)

    for a_domain in a_domains:
        a_domain.set_domain_id(separator_1, separator_2, separator_3)

    return a_domains