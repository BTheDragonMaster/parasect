# -*- coding: utf-8 -*-

"""API for PARASECT."""

import os
from typing import Dict, List, Optional, Union, Any

from sklearn.ensemble import RandomForestClassifier

from parasect.core.constants import FINGERPRINTS_FILE, INCLUDED_SUBSTRATES_FILE, SMILES_FILE, \
    INCLUDED_SUBSTRATES_FILE_BACTERIAL, BACTERIAL_FINGERPRINTS_FILE
from parasect.core.domain import AdenylationDomain
from parasect.database.build_database import AdenylationDomain as ADomain
from parasect.core.featurisation import get_domain_features, get_domains
from parasect.core.parsing import bitvector_from_smiles, data_from_substrate_names, parse_substrate_list

from parasect.core.genbank import fetch_from_genbank
from parasect.core.tabular import Tabular
from parasect.core.parasect_result import Result


def sort_results(results: List["AnnotationResult"]) -> List["ProteinResult"]:
    protein_name_to_results: Dict[str, List[AnnotationResult]] = {}
    protein_name_to_sequence: Dict[str, str] = {}
    for result in results:
        if result.paras_result.domain.protein_name not in protein_name_to_results:
            protein_name_to_results[result.paras_result.domain.protein_name] = []
            protein_name_to_sequence[result.paras_result.domain.protein_name] = result.paras_result.domain.protein_sequence
        protein_name_to_results[result.paras_result.domain.protein_name].append(result)

    protein_results = []
    for protein_name, results in protein_name_to_results.items():
        protein_results.append(ProteinResult(protein_name, protein_name_to_sequence[protein_name], results))

    return protein_results


class ProteinResult:
    """Protein result class."""

    def __init__(
            self,
            protein_name: str,
            sequence: str,
            results: List["AnnotationResult"]
    ) -> None:
        """Initialise the Protein Result class.

        :param protein_name: str, name of the protein
        :param results: List of results for that protein.
        """
        self._protein_name = protein_name
        self._sequence = sequence
        self.results = results

    def to_json(self) -> Dict[str, Union[str, int, List[Dict[str, Union[str, float]]]]]:
        """Return the Result as a JSON serialisable dictionary.

        :return: JSON serialisable dictionary.
        :rtype: Dict[str, Union[str, List[(float, str)]]
        """
        return dict(
            protein_name=self._protein_name,
            sequence=self._sequence,
            results=[result.to_json() for result in self.results]
        )


class AnnotationResult:

    def __init__(self,
                 paras_result: "Result",
                 sequence_matches: List[ADomain],
                 synonym_matches: List[ADomain]):
        self.paras_result = paras_result
        self.sequence_matches: List[dict[str, Any]] = [x.to_json() for x in sequence_matches]
        self.synonym_matches: List[dict[str, Any]] = [y.to_json() for y in synonym_matches]

    def to_json(self):
        return dict(
            paras_result=self.paras_result.to_json(),
            sequence_matches=self.sequence_matches,
            synonym_matches=self.synonym_matches
        )


def run_paras(
    selected_input: str,
    selected_input_type: str,
    path_temp_dir: str,
    model: RandomForestClassifier,
    use_structure_guided_alignment: bool = False,
) -> List[Result]:
    """Predict adenylation domain substrate specificity with PARAS on raw input.

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
    if selected_input_type == "fasta":
        file_name = "input.fasta"
    elif selected_input_type == 'gbk':
        file_name = "input.gbk"
    elif selected_input_type == 'accession':
        file_name = "input.fasta"
    else:
        raise ValueError(f"Unrecognised input type: {selected_input_type}")

    input_file = os.path.join(path_temp_dir, file_name)
    if selected_input_type == 'accession':
        accessions = selected_input.split(';')
        fetch_from_genbank(accessions, input_file)

    else:
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
            smiles = str(smiles_data.get_row_value(row_id=name, column_name="smiles"))
        except KeyError:
            raise KeyError(f"smiles not found for substrate {name}")

        domain_prediction_smiles.append(smiles)

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
    bacterial_only: bool = False
) -> List[Result]:
    """Predict adenylation domain substrate specificity with PARASECT on raw input.

    :param selected_input: Selected input source.
    :type selected_input: str
    :param selected_input_type: Selected input type: fasta or gbk.
    :type selected_input_type: str
    :param path_temp_dir: Path to the temp directory.
    :type path_temp_dir: str
    :param model: Random forest classifier model.
    :type model: RandomForestClassifier
    :param custom_substrate_names: Custom substrate names.
    :type custom_substrate_names: Optional[List[str]]
    :param custom_substrate_smiles: Custom substrate SMILES.
    :type custom_substrate_smiles: Optional[List[str]]
    :param only_custom: Only use custom substrates.
    :type only_custom: bool
    :param use_structure_guided_alignment: Use structure guided alignment.
    :type use_structure_guided_alignment: bool
    :param bacterial_only: Use bacterial model
    :type bacterial_only: bool
    :return: Results.
    :rtype: List[Result]
    :raises ValueError: If substrate names, smiles, and fingerprints are not of the same length.
    :raises RuntimeError: If no feature vectors are found.
    """
    # get names of included substrates
    if bacterial_only:
        substrates_file = INCLUDED_SUBSTRATES_FILE_BACTERIAL
        fingerprints_file = BACTERIAL_FINGERPRINTS_FILE
    else:
        substrates_file = INCLUDED_SUBSTRATES_FILE
        fingerprints_file = FINGERPRINTS_FILE

    included_subs = parse_substrate_list(substrates_file)

    # get smiles and fingerprints for the included substrates
    if not only_custom:

        sub_names, sub_smiles, sub_fps = data_from_substrate_names(included_subs, bacterial_only=bacterial_only)
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
            custom_fp = bitvector_from_smiles(custom_smiles, fingerprints_file)

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

    # parse and return results
    return results


def run_paras_for_signatures(
    domains: List[AdenylationDomain], model: RandomForestClassifier
) -> List[Result]:
    """Predict adenylation domain substrate specificity with PARAS on domains.

    :param domains: List of adenylation domains.
    :type domains: List[AdenylationDomain]
    :param model: Random forest classifier model.
    :type model: RandomForestClassifier
    :return: Results.
    :rtype: List[Result]
    :raises ValueError: If no adenylation domains are supplied.
    :raises KeyError: If smiles not found for substrate.
    """
    # check if list of domains is empty
    if not domains:
        raise ValueError("no adenylation domains supplied")

    # turn signatures into features
    domain_features = [get_domain_features(s.extended_signature) for s in domains]

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
            smiles = str(smiles_data.get_row_value(row_id=name, column_name="smiles"))
        except KeyError:
            raise KeyError(f"smiles not found for substrate {name}")

        domain_prediction_smiles.append(smiles)

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
