# -*- coding: utf-8 -*-

"""API for PARASECT."""

import os
from collections import OrderedDict
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import numpy.typing as npt
from sklearn.ensemble import RandomForestClassifier

from parasect.core.featurisation import domains_to_features, get_domains
from parasect.core.helpers import clear_temp_dir


class Result:
    """Result class."""

    def __init__(
        self,
        domain_id: Optional[str],
        sequence: Optional[str],
        signature: Optional[str],
        extended_signature: Optional[str],
        predictions: Optional[List[Tuple[float, str]]],
    ) -> None:
        """Initialise the Result class.

        :param domain_id: Domain ID.
        :type domain_id: Optional[str]
        :param sequence: Sequence.
        :type sequence: Optional[str]
        :param signature: Signature.
        :type signature: Optional[str]
        :param extended_signature: Extended signature.
        :type extended_signature: Optional[str]
        :param predictions: Predictions.
        :type predictions: Optional[List[Tuple[float, str]]]
        """
        if domain_id is None:
            domain_id = ""

        if sequence is None:
            sequence = ""

        if signature is None:
            signature = ""

        if extended_signature is None:
            extended_signature = ""

        if predictions is None:
            predictions = []

        self.domain_id = domain_id
        self.sequence = sequence
        self.signature = signature
        self.extended_signature = extended_signature
        self.predictions = predictions

    def to_json(self) -> Dict[str, Union[str, Dict[str, Union[str, List[Tuple[float, str]]]]]]:
        """Return the Result as a JSON serialisable dictionary.

        :return: JSON serialisable dictionary.
        :rtype: Dict[str, Union[str, List[(float, str)]]
        """
        return dict(
            domain_id=self.domain_id,
            data=dict(
                sequence=self.sequence,
                signature=self.signature,
                extended_signature=self.extended_signature,
                predictions=self.predictions,
            ),
        )


def get_top_predictions_paras(
    amino_acid_classes: npt.NDArray[np.str],
    probabilities: npt.NDArray[np.float],
    top_n: int,
) -> List[Tuple[float, str]]:
    """Return the top N predictions for a given sequence.

    :param amino_acid_classes: Amino acid classes.
    :type amino_acid_classes: npt.NDArray[np.str]
    :param probabilities: Probabilities.
    :type probabilities: npt.NDArray[np.float]
    :param top_n: Number of predictions to return.
    :type top_n: int
    :return: Top N predictions.
    :rtype: List[Tuple[float, str]]
    """
    # zip probabilities and amino acid classes
    probs_and_aa = []
    for i, probability in enumerate(probabilities):
        probs_and_aa.append((probability, amino_acid_classes[i]))

    # sort by probability.
    probs_and_aa.sort(reverse=True)

    # return top N predictions
    return probs_and_aa[:top_n]


def get_top_predictions_parasect(
    seq_id: str,
    id_to_probabilities: Dict[str, List[Tuple[float, str]]],
    top_n: int,
) -> List[Tuple[float, str]]:
    """Return the top N predictions for a given sequence.

    :param seq_id: Sequence ID.
    :type seq_id: str
    :param id_to_probabilities: ID to probabilities.
    :type id_to_probabilities: Dict[str, List[Tuple[float, str]]]
    :param top_n: Number of predictions to return.
    :type top_n: int
    :return: Top N predictions.
    :rtype: List[Tuple[float, str]]
    """
    probabilities = id_to_probabilities[seq_id]
    probabilities.sort(reverse=True)

    return probabilities[:top_n]


def run_paras(
    selected_input: str,
    selected_input_type: str,
    path_temp_dir: str,
    model: RandomForestClassifier,
    use_structure_guided_alignment: bool = False,
    num_predictions_to_report: int = 3,
    first_separator: str = "|",
    second_separator: str = "_",
    third_separator: str = "-",
) -> List[Result]:
    """Predict adenylation domain substrate specificity with PARAS.

    :param selected_input: Selected input.
    :type selected_input: str
    :param selected_input_type: Selected input type: fasta or gbk.
    :type selected_input_type: str
    :param path_temp_dir: Path to the temp directory.
    :type path_temp_dir: str
    :param model: Random forest classifier model.
    :type model: RandomForestClassifier
    :param use_structure_guided_alignment: Use structure guided alignment.
    :type use_structure_guided_alignment: bool
    :param num_predictions_to_report: Number of predictions to report.
    :type num_predictions_to_report: int
    :param first_separator: First separator.
    :type first_separator: str
    :param second_separator: Second separator.
    :type second_separator: str
    :param third_separator: Third separator.
    :type third_separator: str
    :return: Results.
    :rtype: List[Result]
    :raises Exception: If no feature vectors are found.
    """
    # write selected_input to file in temp folder
    file_name = "input.fasta" if selected_input_type == "fasta" else "input.gbk"
    input_file = os.path.join(path_temp_dir, file_name)
    with open(input_file, "w") as fo:
        fo.write(selected_input)

    # get domains
    a_domains = get_domains(
        path_in=input_file,
        path_temp_dir=path_temp_dir,
        extraction_method="profile" if use_structure_guided_alignment else "hmm",
        file_type=selected_input_type.lower(),
        separator_1=first_separator,
        separator_2=second_separator,
        separator_3=third_separator,
    )
    sequence_ids, feature_vectors = domains_to_features(a_domains)

    # check if feature vectors exist
    if not feature_vectors:
        raise Exception("no feature vectors found")

    # run model and retrieve class predictions
    probabilities = model.predict_proba(feature_vectors)

    results = OrderedDict()
    for sequence_idx, sequence_id in enumerate(sequence_ids):

        probs_and_aa = get_top_predictions_paras(
            amino_acid_classes=model.classes_,
            probabilities=probabilities[sequence_idx],
            top_n=num_predictions_to_report,
        )

        results[sequence_id] = probs_and_aa

    # clean up
    model = None
    clear_temp_dir(path_temp_dir, keep=[".gitkeep"])

    # parse and return results
    return [
        Result(
            domain_id=domain.domain_id,
            sequence=domain.sequence,
            signature=domain.signature,
            extended_signature=domain.extended_signature,
            predictions=results.get(domain.domain_id, None),
        )
        for domain in a_domains
    ]


def run_parasect(
    selected_input: str,
    selected_input_type: str,
    path_temp_dir: str,
    model: RandomForestClassifier,
    substrate_names: List[str],
    substrate_fingerprints: List[List[int]],
    use_structure_guided_alignment: bool = False,
    num_predictions_to_report: int = 3,
    first_separator: str = "|",
    second_separator: str = "_",
    third_separator: str = "-",
) -> List[Result]:
    """Run the PARASECT model.

    :param selected_input: Selected input.
    :type selected_input: str
    :param selected_input_type: Selected input type.
    :type selected_input_type: str
    :param path_temp_dir: Path to the temp directory.
    :type path_temp_dir: str
    :param model: Random forest classifier model.
    :type model: RandomForestClassifier
    :param substrate_names: Substrate names.
    :type substrate_names: List[str]
    :param substrate_fingerprints: Substrate fingerprints.
    :type substrate_fingerprints: List[List[int]]
    :param use_structure_guided_alignment: Use structure guided alignment.
    :type use_structure_guided_alignment: bool
    :param num_predictions_to_report: Number of predictions to report.
    :type num_predictions_to_report: int
    :param first_separator: First separator.
    :type first_separator: str
    :param second_separator: Second separator.
    :type second_separator: str
    :param third_separator: Third separator.
    :type third_separator: str
    :return: Results.
    :rtype: List[Result]
    :raises RuntimeError: If no feature vectors are found.
    :raises ValueError: If interaction values are not 0 or 1.
    """
    # write selected_input to file in temp folder
    file_name = "input.fasta" if selected_input_type == "fasta" else "input.gbk"
    input_file = os.path.join(path_temp_dir, file_name)
    with open(input_file, "w") as fo:
        fo.write(selected_input)

    # get domains
    a_domains = get_domains(
        path_in=input_file,
        path_temp_dir=path_temp_dir,
        extraction_method="profile" if use_structure_guided_alignment else "hmm",
        file_type=selected_input_type.lower(),
        separator_1=first_separator,
        separator_2=second_separator,
        separator_3=third_separator,
    )
    sequence_ids, sequence_feature_vectors = domains_to_features(a_domains)

    # check if feature vectors exist
    if not sequence_feature_vectors:
        raise RuntimeError("no feature vectors found")

    # instantiate results dictionary and set up batch processing
    results = OrderedDict()
    batch_size = 1000
    counter = 0
    start = 0
    end = batch_size

    id_to_probabilities: Dict[str, List[Tuple[float, str]]] = {}
    batch_nr = 1

    # run model and retrieve class predictions in batches
    while start < len(sequence_feature_vectors):

        labels, feature_vectors = [], []
        for i, sequence_feature_vector in enumerate(sequence_feature_vectors[start:end]):
            counter += 1
            for j, fingerprint in enumerate(substrate_fingerprints):
                feature_vector = sequence_feature_vector + fingerprint
                label = (sequence_ids[start + i], substrate_names[j])
                labels.append(label)
                feature_vectors.append(feature_vector)

        start = counter
        end = min([counter + batch_size, len(sequence_feature_vectors)])

        probabilities = model.predict_proba(feature_vectors)
        interaction_labels = model.classes_

        if interaction_labels[0] == 1:
            interaction_index = 0
        elif interaction_labels[1] == 1:
            interaction_index = 1
        else:
            raise ValueError("interaction values must be 0 or 1")

        for i, label in enumerate(labels):
            seq_id, substrate = label
            if seq_id not in id_to_probabilities:
                id_to_probabilities[seq_id] = []

            id_to_probabilities[seq_id].append((probabilities[i][interaction_index], substrate))

        batch_nr += 1

    for sequence_id in sequence_ids:
        results[sequence_id] = get_top_predictions_parasect(
            seq_id=sequence_id,
            id_to_probabilities=id_to_probabilities,
            top_n=num_predictions_to_report,
        )

    # clean up
    model = None
    clear_temp_dir(path_temp_dir, keep=[".gitkeep"])

    # parse and return results
    return [
        Result(
            domain_id=domain.domain_id,
            sequence=domain.sequence,
            signature=domain.signature,
            extended_signature=domain.extended_signature,
            predictions=results.get(domain.domain_id, None),
        )
        for domain in a_domains
    ]
