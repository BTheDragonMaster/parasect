import os
import paras.data.sequence_data.amino_acid_properties
from paras.scripts.parsers.parsers import parse_amino_acid_properties

PROPERTIES_FILE = os.path.join(os.path.dirname(paras.data.sequence_data.amino_acid_properties.__file__),
                               '15_aa_properties_normalised.txt')
PROPERTIES = parse_amino_acid_properties(PROPERTIES_FILE)

ONE_HOT_CATEGORIES = ['A', 'C', 'D', 'E', 'F',
                      'G', 'H', 'I', 'K', 'L',
                      'M', 'N', 'P', 'Q', 'R',
                      'S', 'T', 'V', 'W', 'Y']


def domains_to_features(domains, signature="extended", one_hot=False):
    sequence_ids = []
    feature_vectors = []

    for domain in domains:
        sequence_ids.append(domain.domain_id)

        if signature == "extended":
            if not one_hot:
                feature_vector = get_sequence_features(domain.extended_signature)
            else:
                feature_vector = one_hot_encoding(domain.extended_signature)
        elif signature == "short":
            if not one_hot:
                feature_vector = get_sequence_features(domain.signature)
            else:
                feature_vector = one_hot_encoding(domain.signature)
        else:
            raise ValueError(f"Expected 'extended' or 'short'. Got {signature}.")

        feature_vectors.append(feature_vector)

    return sequence_ids, feature_vectors


def get_sequence_features(sequence):
    features = []

    for aa in sequence:
        properties = PROPERTIES[aa]
        features += properties

    return features


def get_features(id_to_sequence, one_hot=False):
    id_to_features = {}

    for seq_id, sequence in id_to_sequence.items():
        if not one_hot:
            id_to_features[seq_id] = get_sequence_features(sequence)
        else:
            id_to_features[seq_id] = one_hot_encoding(sequence)

    return id_to_features


def one_hot_encoding(sequence):
    feature_vector = []
    for res in sequence:
        assert (res in ONE_HOT_CATEGORIES or res == '-')
        for amino_acid in ONE_HOT_CATEGORIES:

            if res == amino_acid:
                feature_vector.append(1)
            else:
                feature_vector.append(0)
    return feature_vector


def get_sequence_features_bulk(id_to_sequence, one_hot=False):
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


def to_feature_vectors(id_to_features):
    ids = []
    feature_vectors = []
    for seq_id, features in id_to_features.items():
        ids.append(seq_id)
        feature_vectors.append(features)

    return ids, feature_vectors
