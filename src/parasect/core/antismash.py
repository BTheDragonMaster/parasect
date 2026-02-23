from sklearn.ensemble import RandomForestClassifier

from parasect.core.parasect_result import MinimalResult
from parasect.core.featurisation import get_domain_features
from parasect.core.constants import INCLUDED_SUBSTRATES_FILE_PARASECT_BACTERIAL, INCLUDED_SUBSTRATES_FILE_PARASECT
from parasect.core.parsing import data_from_substrate_names, parse_substrate_list


def run_paras_minimal(model: RandomForestClassifier, signatures: list[str], names: list[str]) -> list[MinimalResult]:
    """Run PARAS on signatures

    :param model: parasect model
    :param signatures: 34-amino acid active site signatures
    :param names: domain names corresponding to signatures
    """

    features: list[list[float]] = []
    valid_domain_names: list[str] = []
    valid_signatures: list[str] = []
    for i, signature in enumerate(signatures):
        if len(signature) == 34:
            features.append(get_domain_features(signature))
            valid_domain_names.append(names[i])
            valid_signatures.append(signature)

    prediction_matrix = model.predict_proba(features)

    results: list[MinimalResult] = []

    for i, predictions in enumerate(prediction_matrix):
        result = MinimalResult(valid_domain_names[i], valid_signatures[i], predictions, model.classes_)
        results.append(result)

    return results


def run_parasect_minimal(model: RandomForestClassifier, signatures: list[str], names: list[str],
                         bacterial_only: bool = True) -> list[MinimalResult]:
    """Run PARASECT on signatures

    :param model: parasect model
    :param signatures: 34-amino acid active site signatures
    :param names: domain names corresponding to signatures
    :param bacterial_only: True if bacterial model is used, False otherwise
    """
    if bacterial_only:
        included_substrates_file = INCLUDED_SUBSTRATES_FILE_PARASECT_BACTERIAL
    else:
        included_substrates_file = INCLUDED_SUBSTRATES_FILE_PARASECT

    included_substrates = parse_substrate_list(included_substrates_file)

    substrate_names, substrate_smiles, substrate_fingerprints = data_from_substrate_names(included_substrates,
                                                                                          bacterial_only=bacterial_only)

    valid_domain_features: list[list[float]] = []
    valid_domain_names: list[str] = []
    valid_signatures: list[str] = []

    for i, signature in enumerate(signatures):
        if len(signature) == 34:
            valid_domain_features.append(get_domain_features(signature) )
            valid_domain_names.append(names[i])
            valid_signatures.append(signature)

    results: list[MinimalResult] = []

    for i, domain_features in enumerate(valid_domain_features):
        domain_name = valid_domain_names[i]
        signature = valid_signatures[i]

        full_feature_vectors: list[list[float]] = []

        for j, substrate_fingerprint in enumerate(substrate_fingerprints):

            # concatenate features
            concatenated_features = domain_features + substrate_fingerprint
            full_feature_vectors.append(concatenated_features)

        domain_predictions = model.predict_proba(full_feature_vectors)
        domain_predictions = domain_predictions[:, 1].tolist()

        result = MinimalResult(domain_name, signature, domain_predictions, substrate_names)

        results.append(result)

    return results


