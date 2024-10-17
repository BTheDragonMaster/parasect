# -*- coding: utf-8 -*-

"""Script for training the PARAS model for predicting common substrate specificities."""

import argparse
import logging
import os
import time
from typing import List, Tuple

import numpy as np
from joblib import load, dump
from sklearn.ensemble import RandomForestClassifier

# WOLS870101, WOLS870102, WOLS870103, FAUJ880109, GRAR740102, RADA880108, ZIMJ680103, TSAJ990101, CHOP780201, CHOP780202, CHOP780203, ZIMJ680104, NEU1, NEU2, NEU3
PROPERTIES = {
    "-": [0.00, 0.00, 0.00, 1, 8.3, 0.21, 13.59, 145.2, 1.00, 1.03, 0.99, 6.03, 0.06, 0.00, 0.10],
    "A": [0.07, -1.73, 0.09, 0, 8.1, -0.06, 0.00, 90.0, 1.42, 0.83, 0.66, 6.00, 0.06, -0.25, 0.25],
    "C": [0.71, -0.97, 4.13, 0, 5.5, 1.36, 1.48, 103.3, 0.70, 1.19, 1.19, 5.05, -0.56, -0.40, -0.14],
    "D": [3.64, 1.13, 2.36, 1, 13.0, -0.80, 49.70, 117.3, 1.01, 0.54, 1.46, 2.77, 0.97, -0.08, 0.08],
    "E": [3.08, 0.39, -0.07, 1, 12.3, -0.77, 49.90, 142.2, 1.51, 0.37, 0.74, 3.22, 0.85, -0.10, -0.05],
    "F": [-4.92, 1.30, 0.45, 0, 5.2, 1.27, 0.35, 191.9, 1.13, 1.38, 0.60, 5.48, -0.99, 0.18, 0.15],
    "G": [2.23, -5.36, 0.30, 0, 9.0, -0.41, 0.00, 64.9, 0.57, 0.75, 1.56, 5.97, 0.32, -0.32, 0.28],
    "H": [2.41, 1.74, 1.11, 1, 10.4, 0.49, 51.60, 160.0, 1.00, 0.87, 0.95, 7.59, 0.15, -0.03, -0.10],
    "I": [-4.44, -1.68, -1.03, 0, 5.2, 1.31, 0.13, 163.9, 1.08, 1.60, 0.47, 6.02, -1.00, -0.03, 0.10],
    "K": [2.84, 1.41, -3.14, 2, 11.3, -1.18, 49.50, 167.3, 1.16, 0.74, 1.01, 9.74, 1.00, 0.32, 0.11],
    "L": [-4.19, -1.03, -0.98, 0, 4.9, 1.21, 0.13, 164.0, 1.21, 1.30, 0.59, 5.98, -0.83, 0.05, 0.01],
    "M": [-2.49, -0.27, -0.41, 0, 5.7, 1.27, 1.43, 167.0, 1.45, 1.05, 0.60, 5.74, -0.68, -0.01, 0.04],
    "N": [3.22, 1.45, 0.84, 2, 11.6, -0.48, 3.38, 124.7, 0.67, 0.89, 1.56, 5.41, 0.70, -0.06, 0.17],
    "P": [-1.22, 0.88, 2.23, 0, 8.0, 1.1, 1.58, 122.9, 0.57, 0.55, 1.52, 6.30, 0.45, 0.23, 0.41],
    "Q": [2.18, 0.53, -1.14, 2, 10.5, -0.73, 3.53, 149.4, 1.11, 1.10, 0.98, 5.65, 0.71, -0.02, 0.12],
    "R": [2.88, 2.52, -3.44, 4, 10.5, -0.84, 52.00, 194.0, 0.98, 0.93, 0.95, 10.76, 0.80, 0.19, -0.41],
    "S": [1.96, -1.63, 0.57, 1, 9.2, -0.50, 1.67, 95.4, 0.77, 0.75, 1.43, 5.68, 0.48, -0.15, 0.23],
    "T": [0.92, -2.09, -1.40, 1, 8.6, -0.27, 1.66, 121.5, 0.83, 1.19, 0.96, 5.66, 0.38, -0.10, 0.29],
    "V": [-2.69, -2.53, -1.29, 0, 5.9, 1.09, 0.13, 139.0, 1.06, 1.70, 0.50, 5.96, -0.75, -0.19, 0.03],
    "W": [-4.75, 3.65, 0.85, 1, 5.4, 0.88, 2.10, 228.2, 1.08, 1.37, 0.96, 5.89, -0.57, 0.31, 0.34],
    "Y": [-1.39, 2.32, 0.01, 1, 6.2, 0.33, 1.61, 197.0, 0.69, 1.47, 1.14, 5.66, -0.35, 0.40, -0.02],
}


def get_domain_features(amino_acid_sequence: str) -> List[float]:
    """Return a feature list of NRPSPredictor features from a sequence.

    :param amino_acid_sequence: Amino acid sequence.
    :type amino_acid_sequence: str
    :return: features. List of NRPSPredictor features.
    :rtype: List[float]
    """
    features = []
    for amino_acid_id in amino_acid_sequence:
        properties = PROPERTIES[amino_acid_id]
        features.extend(properties)
    return features


def train_paras_model(path_dataset_in: str,  path_model_out: str) -> None:
    """Train the PARAS model for predicting common substrate specificities.

    :param path_dataset_in: Path to the dataset file.
    :type path_dataset_in: str
    :param path_model_out: Path to save the model.
    :type path_model_out: str
    """
    logger = logging.getLogger(__name__)

    # read dataset file
    signatures, labels = [], []
    with open(path_dataset_in, "r") as fo: 
        fo.readline()  # skip header
        for line in fo:
            _, signature, substrate_specificity = line.strip().split("\t")
            signatures.append(signature)
            labels.append(substrate_specificity)

    # compile the dataset
    X = np.array([get_domain_features(signature) for signature in signatures])
    y = np.array(labels)
    logger.debug(f"compiled dataset with shape: {X.shape}, {y.shape}")

    # instantiate the random forest classifier
    model = RandomForestClassifier(
        n_estimators=1000,
        n_jobs=1,
        oob_score=True,
        random_state=25051989,
        class_weight="balanced",
    )

    # train the model
    logger.debug(f"training the model ...")
    model.fit(X, y)
    logger.debug(f"model training complete with out-of-bag score: {model.oob_score_}")

    # save the model
    logger.debug(f"saving the model ...")
    dump(model, path_model_out)
    logger.debug(f"model saved to: {path_model_out}")


def predict_with_paras(model: RandomForestClassifier, signature: str) -> List[Tuple[str, float]]:
    """Predict the substrate specificity of an A-domain signature using PARAS.

    :param model: PARAS model.
    :type model: RandomForestClassifier
    :param signature: A-domain extended signature (35 amino acids).
    :type signature: str
    :return: Predicted common substrate specificity.
    :rtype: str
    """
    if len(signature) != 34:
        raise ValueError("A-domain extended signature must be 34 amino acids long.")

    features = get_domain_features(signature)
    features = np.array(features).reshape(1, -1)
    predictions = model.predict_proba(features)
    class_labels = model.classes_
    
    # zip the class labels and predictions, and sort by prediction probability
    predictions = list(zip(class_labels, predictions[0]))
    predictions = sorted(predictions, key=lambda x: x[1], reverse=True)

    return predictions


def cli() -> argparse.Namespace:
    """Command line interface for training the PARAS common substrates.
    
    :return: command line arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=str, required=True, help="path to PARAS datset")
    parser.add_argument("--model", type=str, required=True, help="path for output model")
    return parser.parse_args()


def main() -> None:
    """Run the training script."""
    # set up logging
    logging.basicConfig(level="DEBUG")
    logger = logging.getLogger(__name__)
    
    # parse command line arguments
    args = cli()
    logger.debug(f"command line arguments: {args}")

    # save the model
    start_time = time.time()
    train_paras_model(args.dataset, args.model)
    logger.debug(f"training complete in {time.time() - start_time:.2f} seconds")

    # load the model again
    logger.debug("reloading the model ...")
    start_time = time.time()
    model = load(args.model)

    # test the model
    logger.debug(f"testing the model ...")
    signature = "LAKAFDAFVAEGILISAGEVNAYGPTEVTVCATQ"
    expected_top_prediction = "tyrosine"
    predictions = predict_with_paras(model, signature)
    assert predictions[0][0] == expected_top_prediction
    logger.debug(f"model test passed")
    logger.debug(f"loaded model and completed test in {time.time() - start_time:.2f} seconds")

    exit(0)


if __name__ == "__main__":
    main()
