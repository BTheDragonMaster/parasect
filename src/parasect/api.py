# -*- coding: utf-8 -*-

"""API for PARASECT."""

import os
from typing import Dict, List, Optional, Union

from sklearn.ensemble import RandomForestClassifier

from parasect.core.constants import (
    FINGERPRINTS_FILE, 
    INCLUDED_SUBSTRATES_FILE,
    SMILES_FILE,
)
from parasect.core.tabular import Tabular
from parasect.core.domain import AdenylationDomain
from parasect.core.parsing import (
    bitvector_from_smiles, 
    data_from_substrate_names,
    parse_substrate_list,
)
from parasect.core.featurisation import get_domain_features, get_domains
from parasect.core.helpers import clear_temp_dir


class Result:
    """Result class."""

    def __init__(
        self,
        domain: AdenylationDomain,
        predictions: List[float],
        prediction_labels: List[str],
        prediction_smiles: List[str],
    ) -> None:
        """Initialise the Result class.

        :param domain: Adenylation domain.
        :type domain: AdenylationDomain
        :param predictions: Predictions.
        :type predictions: List[float]
        :param prediction_labels: Prediction labels.
        :type prediction_labels: List[str]
        :param prediction_smiles: Prediction SMILES.
        :type prediction_smiles: List[str]
        """
        self._domain = domain
        self._predictions = predictions
        self._prediction_labels = prediction_labels
        self._prediction_smiles = prediction_smiles

    def to_json(self) -> Dict[str, Union[str, int, List[Dict[str, Union[str, float]]]]]:
        """Return the Result as a JSON serialisable dictionary.

        :return: JSON serialisable dictionary.
        :rtype: Dict[str, Union[str, List[(float, str)]]
        """
        return dict(
            domain_name=self._domain.protein_name,
            domain_nr=self._domain.domain_nr,
            domain_start=self._domain.start,
            domain_end=self._domain.end,
            domain_sequence=self._domain.sequence,
            domain_signature=self._domain.signature,
            domain_extended_signature=self._domain.extended_signature,
            predictions=[
                dict(
                    substrate_name=sub_name, 
                    substrate_smiles=sub_smiles, 
                    probability=prob
                )
                for sub_name, sub_smiles, prob in zip(
                    self._prediction_labels, 
                    self._prediction_smiles, 
                    self._predictions,
                )
            ],
        )


def run_paras(
    selected_input: str,
    selected_input_type: str,
    path_temp_dir: str,
    model: RandomForestClassifier,
    use_structure_guided_alignment: bool = False,
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
    :return: Results.
    :rtype: List[Result]
    :raises Exception: If no feature vectors are found.
    :raises KeyError: If smiles not found for substrate.
    """
    # write selected_input to file in temp folder
    file_name = "input.fasta" if selected_input_type == "fasta" else "input.gbk"
    input_file = os.path.join(path_temp_dir, file_name)
    with open(input_file, "w") as fo:
        fo.write(selected_input)

    # get domains
    domains = get_domains(
        path_in=input_file,
        path_temp_dir=path_temp_dir,
        extraction_method="profile" if use_structure_guided_alignment else "hmm",
        file_type=selected_input_type.lower(),
    )

    # check if domains were found
    if not domains:
        raise Exception("no adenylation domains found")

    # get domain features
    domain_features = [get_domain_features(domain.extended_signature) for domain in domains]

    # run model and retrieve class predictions
    domain_predictions = model.predict_proba(domain_features)

    # get prediction classes
    domain_prediction_labels = model.classes_.tolist()

    # find smiles for the labels
    domain_prediction_smiles = []
    smiles_data = Tabular(path_in=SMILES_FILE, separator="\t")
    for name in domain_prediction_labels:

        # get smiles, should be available
        try:
            smiles = smiles_data.get_row_value(row_id=name, column_name="smiles")
        except KeyError:
            raise KeyError(f"smiles not found for substrate {name}")

        domain_prediction_smiles.append(smiles)

    # clean up
    model = None
    clear_temp_dir(path_temp_dir, keep=[".gitkeep"])

    # parse and return results
    return [
        Result(
            domain=domain,
            predictions=domain_predictions[domain_idx].tolist(),
            prediction_labels=domain_prediction_labels,
            prediction_smiles=domain_prediction_smiles,
        )
        for domain_idx, domain in enumerate(domains)
    ]


def run_parasect(
    selected_input: str,
    selected_input_type: str,
    path_temp_dir: str,
    model: RandomForestClassifier,
    custom_substrate_names: Optional[List[str]] = None,
    custom_substrate_smiles: Optional[List[str]] = None,
    only_custom: bool = False,
    use_structure_guided_alignment: bool = False,
) -> List[Result]:
    """Run the PARASECT model.

    :param selected_input: Selected input source.
    :type selected_input: str
    :param selected_input_type: Selected input type: fasta or gbk.
    :type selected_input_type: str
    :param path_temp_dir: Path to the temp directory.
    :type path_temp_dir: str
    :param model: Random forest classifier model.
    :type model: RandomForestClassifier
    :param substrate_names: Substrate names.
    :type substrate_names: List[str]
    :param substrate_smiles: Substrate SMILES.
    :type substrate_smiles: List[str]
    :param only_custom: Only use custom substrates.
    :type only_custom: bool
    :param use_structure_guided_alignment: Use structure guided alignment.
    :type use_structure_guided_alignment: bool
    :return: Results.
    :rtype: List[Result]
    :raises ValueError: If substrate names, smiles, and fingerprints are not of the same length.
    :raises RuntimeError: If no feature vectors are found.
    """
    # get names of included substrates
    included_subs = parse_substrate_list(INCLUDED_SUBSTRATES_FILE)

    # get smiles and fingerprints for the included substrates
    if not only_custom:
        sub_names, sub_smiles, sub_fps = data_from_substrate_names(included_subs)
    else:
        sub_names, sub_smiles, sub_fps = [], [], []

    # if either custom substrate names or smiles are None, set them to empty lists
    if custom_substrate_names is not None and custom_substrate_smiles is not None:

        # check if custom input is same length
        if len(custom_substrate_names) != len(custom_substrate_smiles):
            raise ValueError("substrate names and smiles are not of the same length")

        # loop over custom input and calculate fingerprints
        for custom_name, custom_smiles in zip(custom_substrate_names, custom_substrate_smiles):

            # get substrate fingerprints
            custom_fp = bitvector_from_smiles(custom_smiles, FINGERPRINTS_FILE)

            # add custom data to the included data
            sub_names.append(custom_name)
            sub_smiles.append(custom_smiles)
            sub_fps.append(custom_fp)

    # write selected_input to file in temp folder
    file_name = "input.fasta" if selected_input_type == "fasta" else "input.gbk"
    input_file = os.path.join(path_temp_dir, file_name)
    with open(input_file, "w") as fo:
        fo.write(selected_input)

    # get domains
    domains = get_domains(
        path_in=input_file,
        path_temp_dir=path_temp_dir,
        extraction_method="profile" if use_structure_guided_alignment else "hmm",
        file_type=selected_input_type.lower(),
    )

    # check if domains were found
    if not domains:
        raise RuntimeError("no adenylation domains found")

    # get domain features
    domain_features = [get_domain_features(domain.extended_signature) for domain in domains]

    # compile results
    results = []
    for domain_idx in range(len(domains)):

        # gather all protein+substrate feature vectors
        domain_feature_vectors = []

        # get protein features
        protein_features = domain_features[domain_idx]

        for substrate_idx in range(len(sub_fps)):

            # get substrate fingerprint
            substrate_fingerprint = sub_fps[substrate_idx]

            # concatenate features
            concatenated_features = protein_features + substrate_fingerprint
            domain_feature_vectors.append(concatenated_features)

        # run model and retrieve class predictions
        domain_predictions = model.predict_proba(domain_feature_vectors)
        domain_predictions = domain_predictions[:, 1].tolist()

        # parse and return results
        result = Result(
            domain=domains[domain_idx],
            predictions=domain_predictions,
            prediction_labels=sub_names,
            prediction_smiles=sub_smiles,
        )

        results.append(result)

    # clean up
    model = None
    clear_temp_dir(path_temp_dir, keep=[".gitkeep"])

    # parse and return results
    return results
