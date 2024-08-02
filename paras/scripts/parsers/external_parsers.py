import os
from dataclasses import dataclass

import paras.data.benchmarking
import paras.data.benchmarking.benchmarking_set
import paras.data
from paras.scripts.parsers.tabular import Tabular
from paras.scripts.parsers.parsers import parse_taxonomy
from pprint import pprint

ADENPREDICTOR_PREDICTIONS = os.path.join(os.path.dirname(paras.data.benchmarking.__file__), 'adenpredictor_results.tsv')
ADENPREDICTOR_ABBREVIATIONS = os.path.join(os.path.dirname(paras.data.benchmarking.__file__), 'adenpredictor_abbreviations.txt')

SANDPUMA_IDENTITY = os.path.join(os.path.dirname(paras.data.benchmarking.benchmarking_set.__file__), 'pid.res.tsv')
SANDPUMA_INDIVIDUAL_PREDICTIONS = os.path.join(os.path.dirname(paras.data.benchmarking.benchmarking_set.__file__), 'ind.res.tsv')
SANDPUMA_PREDICTIONS = os.path.join(os.path.dirname(paras.data.benchmarking.benchmarking_set.__file__), 'sandpuma.tsv')
FORCED_PREDICTIONS = os.path.join(os.path.dirname(paras.data.benchmarking.benchmarking_set.__file__), 'ens.res.tsv')

TAXONOMY = os.path.join(os.path.dirname(paras.data.__file__), 'taxonomy.txt')

SANDPUMA_ABBREVIATIONS = os.path.join(os.path.dirname(paras.data.benchmarking.__file__), 'sandpuma_abbreviations.txt')
NRPSPREDICTOR_ABBREVIATIONS = os.path.join(os.path.dirname(paras.data.benchmarking.__file__), 'nrpspredictor_abbreviations.txt')
NRPSPREDICTOR_PREDICTIONS_BACTERIAL = os.path.join(os.path.dirname(paras.data.benchmarking.__file__), 'nrpspredictor_output_bacterial.txt')
NRPSPREDICTOR_PREDICTIONS_FUNGAL = os.path.join(os.path.dirname(paras.data.benchmarking.__file__), 'nrpspredictor_output_fungal.txt')


@dataclass
class AdenPredictorPrediction:
    name: str
    prediction: str
    signature: str


@dataclass
class SandpumaPrediction:
    name: str
    identity: float = None
    predicat_mp: str = None
    predicat_snn: str = None
    asm: str = None
    phmm: str = None
    svm: str = None
    ensemble: str = None
    sandpuma: str = None
    ensemble_score: str = None
    path_accuracy: float = None


@dataclass
class NrpsPredictorPrediction:
    name: str
    prediction: str
    bacterial_prediction: str
    fungal_prediction: str


def parse_identity_file(identity_file):
    name_to_identity = {}
    with open(identity_file, 'r') as identities:
        for line in identities:
            line = line.strip()
            if line:
                domain_name, identity = line.split('\t')
                identity = float(identity)
                name_to_identity[domain_name] = identity

    return name_to_identity


def parse_external_mapping(mapping_file):
    mapping_data = Tabular(mapping_file, [0])
    print(mapping_data.categories)
    print("full name" in mapping_data.categories)
    abbreviation_to_substrate = {}
    for datapoint in mapping_data.data:
        abbreviation = mapping_data.get_value(datapoint, "abbreviation").lower()
        print(abbreviation)
        abbreviation_to_substrate[abbreviation] = \
            mapping_data.get_value(datapoint, "full name")

    return abbreviation_to_substrate


def get_nrpspredictor_substrate(raw_output, abbreviation_to_substrate):
    abbreviation = raw_output.split('(')[0].lower()
    if abbreviation.lower() in abbreviation_to_substrate:
        prediction = abbreviation_to_substrate[abbreviation]
    elif abbreviation.upper() == "N/A":
        prediction = 'N/A'
    else:
        raise Exception(f"Cannot link abbreviation '{abbreviation}' to paras substrate")

    return prediction


def parse_adenpredictor_predictions(prediction_file, abbreviation_to_substrate):
    name_to_prediction = {}
    adenpredictor_data = Tabular(prediction_file, [0])
    for datapoint in adenpredictor_data.data:
        domain_name = adenpredictor_data.get_value(datapoint, 'Id')
        signature = adenpredictor_data.get_value(datapoint, 'Signature')
        prediction_abbreviation = adenpredictor_data.get_value(datapoint, 'Prediction 1')
        prediction = abbreviation_to_substrate[prediction_abbreviation]
        name_to_prediction[domain_name] = AdenPredictorPrediction(domain_name, prediction, signature)

    return name_to_prediction


def parse_nrpspredictor_predictions(prediction_file_bacterial, prediction_file_fungal, abbreviation_to_substrate):
    name_to_prediction = {}
    domain_to_taxonomy = parse_taxonomy(TAXONOMY)
    nrpspredictor_data_bacterial = Tabular(prediction_file_bacterial, [0])
    nrpspredictor_data_fungal = Tabular(prediction_file_fungal, [0])
    for datapoint in nrpspredictor_data_bacterial.data:

        domain_name = nrpspredictor_data_bacterial.get_value(datapoint, "Name")
        abbreviated_fungal_prediction = nrpspredictor_data_fungal.get_value(datapoint, "SingleV2")
        abbreviated_bacterial_prediction = nrpspredictor_data_bacterial.get_value(datapoint, "SingleV2")
        fungal_prediction = get_nrpspredictor_substrate(abbreviated_fungal_prediction, abbreviation_to_substrate)
        bacterial_prediction = get_nrpspredictor_substrate(abbreviated_bacterial_prediction, abbreviation_to_substrate)

        if domain_name not in domain_to_taxonomy:
            prediction = bacterial_prediction
        elif 'Fungi' in domain_to_taxonomy[domain_name]:
            prediction = fungal_prediction
        else:
            prediction = bacterial_prediction

        name_to_prediction[domain_name] = NrpsPredictorPrediction(domain_name, prediction, bacterial_prediction,
                                                                  fungal_prediction)
    return name_to_prediction


def parse_prediction(prediction_file, abbreviation_to_substrate):
    name_to_method_to_prediction = {}
    with open(prediction_file, 'r') as predictions:
        for line in predictions:
            line_data = line.split('\t')

            domain_name, method, result = line_data[:3]
            result = result.strip()
            if domain_name not in name_to_method_to_prediction:
                name_to_method_to_prediction[domain_name] = {}
            if len(line_data) >= 4:
                ensemble_prediction_score = line_data[3]
                if ensemble_prediction_score.strip():
                    if '=' not in ensemble_prediction_score:
                        name_to_method_to_prediction[domain_name]['score'] = float(ensemble_prediction_score.strip())
                    else:
                        ensemble_prediction_score = ensemble_prediction_score.strip().split('=')[1]
                        name_to_method_to_prediction[domain_name]['score'] = float(ensemble_prediction_score)

            if len(line_data) == 5:
                path_accuracy = line_data[4].split('=')[1].strip()
                if path_accuracy:
                    name_to_method_to_prediction[domain_name]['path_score'] = float(path_accuracy)
                else:
                    name_to_method_to_prediction[domain_name]['path_score'] = None

            results = result.split('|')
            for result in results:
                name_to_method_to_prediction[domain_name][method] = []

                if result.lower() in abbreviation_to_substrate:

                    name_to_method_to_prediction[domain_name][method].append(abbreviation_to_substrate[result.lower()])
                elif result in {"no_confident_result", "no_force_needed", "no_call", "N/A"}:
                    name_to_method_to_prediction[domain_name][method].append(result)
                else:
                    raise Exception(f"Cannot link abbreviation '{result}' to paras substrate")
    return name_to_method_to_prediction


def parse_nrpspredictor_data():
    abbreviation_to_substrate = parse_external_mapping(NRPSPREDICTOR_ABBREVIATIONS)
    domain_to_prediction = parse_nrpspredictor_predictions(NRPSPREDICTOR_PREDICTIONS_BACTERIAL,
                                                           NRPSPREDICTOR_PREDICTIONS_FUNGAL, abbreviation_to_substrate)
    return domain_to_prediction


def parse_adenpredictor_data():
    abbreviation_to_substrate = parse_external_mapping(ADENPREDICTOR_ABBREVIATIONS)
    domain_to_prediction = parse_adenpredictor_predictions(ADENPREDICTOR_PREDICTIONS, abbreviation_to_substrate)

    return domain_to_prediction


def parse_sandpuma_data():
    abbreviation_to_substrate = parse_external_mapping(SANDPUMA_ABBREVIATIONS)
    individual_predictions = parse_prediction(SANDPUMA_INDIVIDUAL_PREDICTIONS, abbreviation_to_substrate)
    ensemble_predictions = parse_prediction(FORCED_PREDICTIONS, abbreviation_to_substrate)
    sandpuma_predictions = parse_prediction(SANDPUMA_PREDICTIONS, abbreviation_to_substrate)

    domain_to_id = parse_identity_file(SANDPUMA_IDENTITY)
    domain_to_prediction = {}

    for domain, identity in domain_to_id.items():
        sandpuma_prediction = SandpumaPrediction(domain, identity=identity)
        sandpuma_prediction.sandpuma = sandpuma_predictions[domain]["SANDPUMA"]
        sandpuma_prediction.ensemble = ensemble_predictions[domain]["ENS"]
        sandpuma_prediction.ensemble_score = sandpuma_predictions[domain]["score"]
        sandpuma_prediction.path_accuracy = sandpuma_predictions[domain]["path_score"]
        sandpuma_prediction.asm = individual_predictions[domain]["ASM"]
        sandpuma_prediction.svm = individual_predictions[domain]["SVM"]
        sandpuma_prediction.phmm = individual_predictions[domain]["pHMM"]
        sandpuma_prediction.predicat_mp = individual_predictions[domain]["prediCAT_MP"]
        sandpuma_prediction.predicat_snn = individual_predictions[domain]["prediCAT_SNN"]

        domain_to_prediction[domain] = sandpuma_prediction

    return domain_to_prediction


if __name__ == "__main__":
    sandpuma_domain_to_prediction = parse_sandpuma_data()
    pprint(sandpuma_domain_to_prediction)

    nrpspredictor_domain_to_prediction = parse_nrpspredictor_data()
    pprint(nrpspredictor_domain_to_prediction)

    adenpredictor_domain_to_prediction = parse_adenpredictor_data()
    pprint(adenpredictor_domain_to_prediction)


