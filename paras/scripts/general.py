import os
from typing import Optional, List
from collections import OrderedDict

from joblib import load

import paras.data.compound_data

from paras.scripts.parsers.fasta import yield_fasta
from paras.scripts.parsers.gbk_to_fasta import proteins_from_genbank
from paras.scripts.feature_extraction.sequence_feature_extraction.extract_domains import domains_from_fasta
from paras.scripts.feature_extraction.sequence_feature_extraction.sequence_labels import SEPARATOR_1, \
    SEPARATOR_2, SEPARATOR_3
from paras.scripts.feature_extraction.sequence_feature_extraction.seq_to_features import domains_to_features, \
    get_sequence_features, to_feature_vectors
from paras.scripts.feature_extraction.sequence_feature_extraction.rename_sequences import rename_sequences, \
    reverse_renaming
from paras.scripts.data_processing.temp import clear_temp
from paras.scripts.parsers.parsers import parse_substrate_list
from paras.scripts.feature_extraction.compound_feature_extraction.fingerprinting import bitvector_from_smiles, \
    bitvectors_from_substrate_names
from paras.scripts.feature_extraction.sequence_feature_extraction.adenylation_domain import VALID_CHARACTERS

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


def parasect_from_extended_signature(extended_signature, verbose=False):
    if not all([char in VALID_CHARACTERS for char in extended_signature]):
        return None
    else:
        pass


def paras_from_extended_signatures(extended_signatures, model_path, nr_results=3, verbose=False):
    feature_vectors = []
    valid_indices = []

    if verbose:
        print("Loading PARAS classifier..")

    classifier = load(model_path) # PARAS model

    for i, extended_signature in enumerate(extended_signatures):
        if all([char in VALID_CHARACTERS for char in extended_signature]):
            if verbose:
                print("Extracting features..")

            feature_vector = get_sequence_features(extended_signature)
            feature_vectors.append(feature_vector)
            valid_indices.append(i)

    if verbose:
        print("Predicting substrate..")
    probabilities = classifier.predict_proba(feature_vectors)
    amino_acid_classes = classifier.classes_

    result_index = 0
    predictions = []

    for i in range(len(extended_signatures)):
        if i in valid_indices:
            probability_list = probabilities[result_index]
            probs_and_aa = get_top_n_aa_paras(amino_acid_classes, probability_list, nr_results)
            predictions.append(probs_and_aa)
            result_index += 1
        else:
            predictions.append(None)

    return predictions


def get_domains(input_file, extraction_method, job_name, separator_1, separator_2, separator_3, verbose, file_type, temp_dir):
    assert extraction_method in ['hmm', 'profile'], f"Only supported extraction methods are 'hmm' or 'profile'. Got {extraction_method}."
    assert file_type in ['fasta', 'gbk'], f"Only supported file types are 'fasta' or 'gbk'. Got {file_type}."
    
    if file_type == 'gbk':
        original_fasta = os.path.join(temp_dir, 'proteins_from_genbank.fasta')
        proteins_from_genbank(input_file, original_fasta) # TODO: ???
        mapping_file, renamed_fasta_file = rename_sequences(original_fasta, temp_dir)
    else:
        mapping_file, renamed_fasta_file = rename_sequences(input_file, temp_dir)

    if extraction_method == 'profile':
        a_domains = domains_from_fasta(renamed_fasta_file, temp_dir=temp_dir, job_name=job_name, profile=True, verbose=verbose)
    elif extraction_method == 'hmm':
        a_domains = domains_from_fasta(renamed_fasta_file, temp_dir=temp_dir, job_name=job_name, verbose=verbose)
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


def run_paras(
    fasta_file, 
    job_name, 
    model_path_paras,
    model_path_paras_all,
    model_path_paras_onehot,
    temp_dir,
    hmmer_path,
    extraction_method='hmm', 
    save_signatures=False,
    save_extended_signatures=False, 
    save_domains=False, 
    nr_results=3, 
    separator_1=SEPARATOR_1,
    separator_2=SEPARATOR_2, 
    separator_3=SEPARATOR_3, 
    verbose=False, 
    file_type='fasta', 
    one_hot=False,
    all_substrates=False,
):
    assert file_type in ['fasta', 'gbk']

    a_domains = get_domains(fasta_file, extraction_method, job_name, separator_1, separator_2, separator_3, verbose,
                            file_type, temp_dir, hmmer_path)

    if verbose:
        print("Extracting features..")
    ids, feature_vectors = domains_to_features(a_domains, one_hot=one_hot)

    results = OrderedDict()
    if feature_vectors:
        if verbose:
            print("Loading PARAS classifier..")
        if not one_hot:
            if not all_substrates:
                assert model_path_paras is not None
                classifier = load(model_path_paras)
            else:
                assert model_path_paras_all is not None
                classifier = load(model_path_paras_all)
        else:
            assert model_path_paras_onehot is not None
            classifier = load(model_path_paras_onehot)
        if verbose:
            print("Predicting substrates..")
        probabilities = classifier.predict_proba(feature_vectors)
        amino_acid_classes = classifier.classes_

        for i, seq_id in enumerate(ids):
            probability_list = probabilities[i]
            probs_and_aa = get_top_n_aa_paras(amino_acid_classes, probability_list, nr_results)
            results[seq_id] = probs_and_aa

    else:
        if verbose:
            print("No A domains found.")
        clear_temp(temp_dir)
        return results

    if verbose:
        print("Writing results..")  

    clear_temp(temp_dir)

    domain_results = {}
    for domain in a_domains:
        domain_results[domain.domain_id] = {}
        if save_domains:
            domain_results[domain.domain_id]['sequence'] = domain.sequence
        if save_signatures:
            domain_results[domain.domain_id]['signature'] = domain.signature
        if save_extended_signatures:
            domain_results[domain.domain_id]['extended_signature'] = domain.extended_signature

    for domain_id in results:
        preds = results[domain_id]
        domain_results[domain_id]['predictions'] = preds
        
    return domain_results

def run_parasect(
    fasta_file, 
    job_name, 
    model_path_parasect,
    model_path_parasect_onehot,
    temp_dir,
    hmmer_path,
    exclude_default_substrates=False, 
    smiles: Optional[List] = None,
    substrate_names: Optional[List] = None, 
    extraction_method='hmm', 
    save_signatures=False,
    save_extended_signatures=False, 
    save_domains=False, 
    nr_results=3, 
    separator_1=SEPARATOR_1,
    separator_2=SEPARATOR_2, 
    separator_3=SEPARATOR_3, 
    verbose=False, 
    file_type='fasta', 
    one_hot=False
):

    assert file_type in ['fasta', 'gbk']

    a_domains = get_domains(fasta_file, extraction_method, job_name, separator_1, separator_2, separator_3, verbose,
                            file_type, temp_dir, hmmer_path)

    if verbose:
        print("Extracting domain features..")
    sequence_ids, sequence_feature_vectors = domains_to_features(a_domains, one_hot=one_hot)

    if verbose:
        print("Extracting substrate features..")

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

    results = OrderedDict()
    if fingerprints and sequence_feature_vectors:
        if verbose:
            print("Loading PARASECT classifier..")

        if not one_hot:
            classifier = load(model_path_parasect)
        else:
            classifier = load(model_path_parasect_onehot)

        if verbose:
            print("Predicting substrates..")

        batch_size = 1000

        counter = 0
        start = 0
        end = batch_size

        id_to_probabilities = {}

        batch_nr = 1

        while start < len(sequence_feature_vectors):

            if verbose:
                print(f"\tPredicting batch {batch_nr}..")

            labels, feature_vectors = [], []
            for i, sequence_feature_vector in enumerate(sequence_feature_vectors[start:end]):
                counter += 1
                for j, fingerprint in enumerate(fingerprints):
                    feature_vector = sequence_feature_vector + fingerprint
                    label = (sequence_ids[start + i], substrates[j])
                    labels.append(label)
                    feature_vectors.append(feature_vector)

            start = counter
            end = min([counter + batch_size, len(sequence_feature_vectors)])

            probabilities = classifier.predict_proba(feature_vectors)
            interaction_labels = classifier.classes_

            if interaction_labels[0] == 1:
                interaction_index = 0
            elif interaction_labels[1] == 1:
                interaction_index = 1
            else:
                raise ValueError("Interaction values must be 0 and 1")

            for i, label in enumerate(labels):
                seq_id, substrate = label
                if seq_id not in id_to_probabilities:
                    id_to_probabilities[seq_id] = []

                id_to_probabilities[seq_id].append((probabilities[i][interaction_index], substrate))

            batch_nr += 1

        for seq_id in sequence_ids:
            results[seq_id] = get_top_n_aa_parasect(seq_id, id_to_probabilities, nr_results)

    else:
        if verbose:
            print("No A domains found.")
        clear_temp()
        return results

    if verbose:
        print("Writing results..")

    # write_files(a_domains, out_dir, job_name, save_domains, save_signatures, save_extended_signatures)
    clear_temp(temp_dir)

    return results


def run_custom_model(fasta_file: str, model_path: str, nr_results: int = 3, seq_len: int = 33, verbose: bool = False):
    """
    Returns a dictionary of paras results

    Parameters
    ----------
    fasta_file: str, path to fasta file. All sequences must have equal length
    model_path: str, path to scikit learn model to run
    nr_results: int, number of top predictions to return
    seq_len: int, length of the sequences in the fasta file
    verbose: bool, print progress if True. Default: False

    Returns
    -------
    results: dict of {seq_id: [(prediction_probability, prediction), ->], ->}, with seq_id str, prediction_probability
        float, and prediction str
    """

    if verbose:
        print("Loading data..")
    id_to_features = {}

    for domain, seq in yield_fasta(fasta_file):
        new_seq = seq.replace('X', '-')
        if all([char in VALID_CHARACTERS for char in new_seq]):
            if len(new_seq) == seq_len:
                id_to_features[domain] = get_sequence_features(new_seq)

            else:
                print(f"Cannot make prediction for domain {domain}. Sequence {seq} is not of length {seq_len}.")
        else:
            print(f"Cannot make prediction for domain {domain}. Sequence {seq} contains invalid characters.")

    ids, feature_vectors = to_feature_vectors(id_to_features)

    results = OrderedDict()
    if feature_vectors:
        if verbose:
            print("Loading classifier..")
        classifier = load(model_path)
        if verbose:
            print("Predicting substrates..")

        probabilities = classifier.predict_proba(feature_vectors)
        amino_acid_classes = classifier.classes_

        for i, seq_id in enumerate(ids):
            probability_list = probabilities[i]
            probs_and_aa = get_top_n_aa_paras(amino_acid_classes, probability_list, nr_results)
            results[seq_id] = probs_and_aa

    else:
        if verbose:
            print("No A domains found.")
        clear_temp()
        return results

    return results
