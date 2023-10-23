import os
from typing import Optional, List
from collections import OrderedDict

from joblib import load

import paras.models.random_forest
import paras.data.compound_data

from paras.scripts.feature_extraction.sequence_feature_extraction.extract_domains import domains_from_fasta
from paras.scripts.feature_extraction.sequence_feature_extraction.sequence_labels import SEPARATOR_1, \
    SEPARATOR_2, SEPARATOR_3
from paras.scripts.feature_extraction.sequence_feature_extraction.seq_to_features import domains_to_features
from paras.scripts.feature_extraction.sequence_feature_extraction.rename_sequences import rename_sequences, \
    reverse_renaming
from paras.scripts.data_processing.temp import TEMP_DIR, clear_temp
from paras.scripts.parsers.parsers import parse_substrate_list
from paras.scripts.feature_extraction.compound_feature_extraction.fingerprinting import bitvector_from_smiles, \
    bitvectors_from_substrate_names


PARAS = os.path.join(os.path.dirname(paras.models.random_forest.__file__), 'class_sequence.paras')
PARASECT = os.path.join(os.path.dirname(paras.models.random_forest.__file__), 'class_sequence.parasect')
INCLUDED_SUBSTRATES = os.path.join(os.path.dirname(paras.data.compound_data.__file__), 'included_substrates.txt')
FINGERPRINTS = os.path.join(os.path.dirname(paras.data.compound_data.__file__), 'fingerprints.txt')


def write_results(results, out_file):
    with open(out_file, 'w') as out:
        out.write(f"sequence_id")
        header_written = False
        for seq_id, probabilities_and_substrates in results.items():
            if not header_written:
                for i in range(len(probabilities_and_substrates)):
                    out.write(f'\tsubstrate_{i + 1}\tconfidence_score_{i + 1}')
                out.write('\n')
            header_written = True

            out.write(f"{seq_id}")
            for probability, substrate in probabilities_and_substrates:
                out.write(f"\t{substrate}\t{probability:.2f}")
            out.write('\n')


def get_top_n_aa_paras(amino_acid_classes, probabilities, n):
    probs_and_aa = []
    for i, probability in enumerate(probabilities):
        probs_and_aa.append((probability, amino_acid_classes[i]))

    probs_and_aa.sort(reverse=True)

    return probs_and_aa[:n]


def get_top_n_aa_parasect(seq_id, id_to_probabilities, n):
    probabilities = id_to_probabilities[seq_id]
    probabilities.sort(reverse=True)
    return probabilities[:n]


def get_domains(fasta_file, extraction_method, job_name, separator_1, separator_2, separator_3):
    assert extraction_method in ['hmm', 'profile']

    mapping_file, renamed_fasta_file = rename_sequences(fasta_file, TEMP_DIR)

    if extraction_method == 'profile':
        a_domains = domains_from_fasta(renamed_fasta_file, job_name=job_name, profile=True)
    elif extraction_method == 'hmm':
        a_domains = domains_from_fasta(renamed_fasta_file, job_name=job_name)
    else:
        raise ValueError(f"Only supported extraction methods are 'hmm' or 'profile'. Got {extraction_method}.")

    reverse_renaming(a_domains, mapping_file)

    for a_domain in a_domains:
        a_domain.set_domain_id(separator_1, separator_2, separator_3)

    return a_domains


def write_files(a_domains, out_dir, job_name, save_domains, save_signatures, save_extended_signatures):
    if save_domains:
        fasta_out = os.path.join(out_dir, f"{job_name}_domains.fasta")
        with open(fasta_out, 'w') as fasta:
            for domain in a_domains:
                domain.write_sequence(fasta)
    if save_signatures:
        fasta_out = os.path.join(out_dir, f"{job_name}_signatures.fasta")
        with open(fasta_out, 'w') as fasta:
            for domain in a_domains:
                domain.write_sequence(fasta, sequence_type='signature')
    if save_extended_signatures:
        fasta_out = os.path.join(out_dir, f"{job_name}_extended_signatures.fasta")
        with open(fasta_out, 'w') as fasta:
            for domain in a_domains:
                domain.write_sequence(fasta, sequence_type='extended_signature')


def run_paras(fasta_file, out_dir, job_name, extraction_method='hmm', save_signatures=False,
              save_extended_signatures=False, save_domains=False, nr_results=3, separator_1=SEPARATOR_1,
              separator_2=SEPARATOR_2, separator_3=SEPARATOR_3):

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    a_domains = get_domains(fasta_file, extraction_method, job_name, separator_1, separator_2, separator_3)
    ids, feature_vectors = domains_to_features(a_domains)

    results = OrderedDict()
    if feature_vectors:
        classifier = load(PARAS)
        probabilities = classifier.predict_proba(feature_vectors)
        amino_acid_classes = classifier.classes_

        for i, seq_id in enumerate(ids):
            probability_list = probabilities[i]
            probs_and_aa = get_top_n_aa_paras(amino_acid_classes, probability_list, nr_results)
            results[seq_id] = probs_and_aa

    else:
        print("No A domains found.")

    write_files(a_domains, out_dir, job_name, save_domains, save_signatures, save_extended_signatures)
    clear_temp()

    return results


def run_parasect(fasta_file, out_dir, job_name, exclude_default_substrates=False, smiles: Optional[List] = None,
                 substrate_names: Optional[List] = None, extraction_method='hmm', save_signatures=False,
                 save_extended_signatures=False, save_domains=False, nr_results=3, separator_1=SEPARATOR_1,
                 separator_2=SEPARATOR_2, separator_3=SEPARATOR_3):

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    a_domains = get_domains(fasta_file, extraction_method, job_name, separator_1, separator_2, separator_3)
    sequence_ids, sequence_feature_vectors = domains_to_features(a_domains)

    if not exclude_default_substrates:
        substrates = parse_substrate_list(INCLUDED_SUBSTRATES)
    else:
        substrates = []

    if substrate_names:
        substrates += substrate_names

    substrates, fingerprints = bitvectors_from_substrate_names(substrates, FINGERPRINTS)

    if smiles:
        for smiles_string in smiles:
            fingerprint = bitvector_from_smiles(smiles_string, FINGERPRINTS)
            fingerprints.append(fingerprint)
            substrates.append(smiles_string)

    labels, feature_vectors = [], []

    for i, sequence_feature_vector in enumerate(sequence_feature_vectors):
        for j, fingerprint in enumerate(fingerprints):
            feature_vector = sequence_feature_vector + fingerprint
            label = (sequence_ids[i], substrates[j])
            labels.append(label)
            feature_vectors.append(feature_vector)

    results = OrderedDict()
    if feature_vectors:
        classifier = load(PARASECT)
        probabilities = classifier.predict_proba(feature_vectors)
        interaction_labels = classifier.classes_

        if interaction_labels[0] == 1:
            interaction_index = 0
        elif interaction_labels[1] == 1:
            interaction_index = 1
        else:
            raise ValueError("Interaction values must be 0 and 1")

        id_to_probabilities = {}

        for i, label in enumerate(labels):
            seq_id, substrate = label
            if seq_id not in id_to_probabilities:
                id_to_probabilities[seq_id] = []

            id_to_probabilities[seq_id].append((probabilities[i][interaction_index], substrate))

        for seq_id in sequence_ids:
            results[seq_id] = get_top_n_aa_parasect(seq_id, id_to_probabilities, nr_results)

    else:
        print("No A domains/substrates found.")

    write_files(a_domains, out_dir, job_name, save_domains, save_signatures, save_extended_signatures)
    clear_temp()

    return results
