import os
from dataclasses import dataclass

import paras.data.benchmarking
import paras.data.benchmarking.benchmarking_set
from paras.scripts.parsers.tabular import Tabular
from pprint import pprint


SANDPUMA_IDENTITY = os.path.join(os.path.dirname(paras.data.benchmarking.benchmarking_set.__file__), 'pid.res.tsv')
SANDPUMA_INDIVIDUAL_PREDICTIONS = os.path.join(os.path.dirname(paras.data.benchmarking.benchmarking_set.__file__), 'ind.res.tsv')
SANDPUMA_PREDICTIONS = os.path.join(os.path.dirname(paras.data.benchmarking.benchmarking_set.__file__), 'sandpuma.tsv')
FORCED_PREDICTIONS = os.path.join(os.path.dirname(paras.data.benchmarking.benchmarking_set.__file__), 'ens.res.tsv')

SANDPUMA_ABBREVIATIONS = os.path.join(os.path.dirname(paras.data.benchmarking.__file__), 'sandpuma_abbreviations.txt')


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


def parse_sandpuma_mapping(sandpuma_mapping_file):
    sandpuma_mapping_data = Tabular(sandpuma_mapping_file, [0])
    abbreviation_to_substrate = {}
    for datapoint in sandpuma_mapping_data.data:
        abbreviation_to_substrate[sandpuma_mapping_data.get_value(datapoint, "abbreviation").lower()] = \
            sandpuma_mapping_data.get_value(datapoint, "full name")

    return abbreviation_to_substrate


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


def parse_sandpuma_data():
    abbreviation_to_substrate = parse_sandpuma_mapping(SANDPUMA_ABBREVIATIONS)
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
    domain_to_prediction = parse_sandpuma_data()
    pprint(domain_to_prediction)


