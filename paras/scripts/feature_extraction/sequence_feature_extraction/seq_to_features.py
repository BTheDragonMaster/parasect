import os
from typing import TYPE_CHECKING, Union
import paras.data.sequence_data.amino_acid_properties
from paras.scripts.parsers.parsers import parse_amino_acid_properties

if TYPE_CHECKING:
    from paras.scripts.feature_extraction.sequence_feature_extraction.adenylation_domain import AdenylationDomain

PROPERTIES_FILE: str = os.path.join(os.path.dirname(paras.data.sequence_data.amino_acid_properties.__file__),
                                    '15_aa_properties_normalised.txt')
PROPERTIES: dict[str, list[float]] = parse_amino_acid_properties(PROPERTIES_FILE)

ONE_HOT_CATEGORIES = ['A', 'C', 'D', 'E', 'F',
                      'G', 'H', 'I', 'K', 'L',
                      'M', 'N', 'P', 'Q', 'R',
                      'S', 'T', 'V', 'W', 'Y']


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


def get_features(id_to_sequence: dict[str, str], one_hot: bool = False) -> dict[str, Union[list[float], list[int]]]:
    """
    Return a dictionary of sequence ids pointing to corresponding feature vectors

    Parameters
    ----------
    id_to_sequence: dict of {sequence_id: sequence, ->}, with sequence_id and sequence both str. Sequence must be an
        amino acid sequence
    one_hot: bool, use one-hot encoding if True, the 15 NRPSPredictor features otherwise

    Returns
    -------
    id_to_features: dict of {sequence_id: [feature, ->], ->}, with sequence_id str and feature int or float

    """
    id_to_features: dict[str, Union[list[int], list[float]]] = {}

    for seq_id, sequence in id_to_sequence.items():
        if not one_hot:
            id_to_features[seq_id] = get_sequence_features(sequence)
        else:
            id_to_features[seq_id] = one_hot_encoding(sequence)

    return id_to_features


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


def get_sequence_features_bulk(id_to_sequence: dict[str, str], one_hot: bool = False):
    """
    Return dictionary of sequence id to feature vector and a list of categories mapping to each feature in the vector

    Parameters
    ----------
    id_to_sequence: dict of {sequence_id: sequence, ->}, with sequence_id and sequence both str. Sequence must be an
        amino acid sequence
    one_hot: bool, use one-hot encoding if True, the 15 NRPSPredictor features otherwise

    Returns
    -------
    id_to_features: dict of {sequence_id: [feature, ->], ->}, with sequence_id str and feature int or float
    all_categories: list of [category, ->], with category string. List has the same length as each feature vector

    """
    id_to_features = {}
    all_categories = []
    seq_length = 0

    if not one_hot:
        _, categories = parse_amino_acid_properties(PROPERTIES_FILE, return_categories=True)

    else:
        categories = ONE_HOT_CATEGORIES

    for seq_id, sequence in id_to_sequence.items():
        seq_length = len(sequence)
        if not one_hot:
            id_to_features[seq_id] = get_sequence_features(sequence)
        else:
            id_to_features[seq_id] = one_hot_encoding(sequence)

    for i in range(1, seq_length + 1):
        for category in categories:
            all_categories.append(f"{category}_res_{i}")

    return id_to_features, all_categories


def to_feature_vectors(id_to_features: dict[str, Union[list[float], list[int]]]) \
        -> tuple[list[str], list[Union[list[float], list[int]]]]:
    """
    Unzip dictionary of sequence ids and corresponding feature vectors into lists of sequence ids and feature vectors

    Parameters
    ----------
    id_to_features: dict of {sequence_id: [feature, ->], ->}, with sequence_id str and feature float or int

    Returns
    -------
    ids: list of [sequence_id, ->], with sequence_id str
    feature_vectors: list of [[feature, ->], ->], with feature int or float. Same length as ids.

    """
    ids: list[str] = []
    feature_vectors: list[Union[list[float], list[int]]] = []
    for seq_id, features in id_to_features.items():
        ids.append(seq_id)
        feature_vectors.append(features)

    return ids, feature_vectors
