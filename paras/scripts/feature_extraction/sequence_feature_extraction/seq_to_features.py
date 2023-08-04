import os
import paras.data.sequence_data.amino_acid_properties
from paras.scripts.parsers.parsers import parse_amino_acid_properties

PROPERTIES_FILE = os.path.join(os.path.dirname(paras.data.sequence_data.amino_acid_properties.__file__), '15_aa_properties_normalised.txt')
PROPERTIES = parse_amino_acid_properties(PROPERTIES_FILE)


def get_sequence_features(sequence):
    features = []

    for aa in sequence:
        properties = PROPERTIES[aa]
        features += properties

    return features


def get_features(id_to_sequence):
    id_to_features = {}

    for seq_id, sequence in id_to_sequence.items():
        id_to_features[seq_id] = get_sequence_features(sequence)

    return id_to_features


def get_sequence_features_bulk(id_to_sequence):
    _, categories = parse_amino_acid_properties(PROPERTIES_FILE, return_categories=True)
    id_to_features = {}
    seq_length = 0

    for seq_id, sequence in id_to_sequence.items():
        seq_length = len(sequence)
        id_to_features[seq_id] = get_sequence_features(sequence)

    all_categories = []

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
