import os
from joblib import load
import paras.models.random_forest
import paras.data.compound_data
from paras.scripts.feature_extraction.sequence_feature_extraction.extract_domains import domains_from_fasta
from paras.scripts.feature_extraction.sequence_feature_extraction.seq_to_features import domains_to_features
from paras.scripts.feature_extraction.sequence_feature_extraction.rename_sequences import rename_sequences, \
    reverse_renaming
from paras.scripts.data_processing.temp import TEMP_DIR, clear_temp
from paras.scripts.feature_extraction.sequence_feature_extraction.sequence_labels import SEPARATOR_1, \
    SEPARATOR_2, SEPARATOR_3

MODEL = os.path.join(os.path.dirname(paras.models.random_forest.__file__), 'class_sequence.paras')


def get_top_x_aa(amino_acid_classes, probabilities, x):
    probs_and_aa = []
    for i, probability in enumerate(probabilities):
        probs_and_aa.append((probability, amino_acid_classes[i]))

    probs_and_aa.sort(reverse=True)

    return probs_and_aa[:x]


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


def run_paras(fasta_file, out_dir, job_name, extraction_method='hmm', save_signatures=False,
              save_extended_signatures=False, save_domains=False, nr_results=3, separator_1=SEPARATOR_1,
              separator_2=SEPARATOR_2, separator_3=SEPARATOR_3):

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    a_domains = get_domains(fasta_file, extraction_method, job_name, separator_1, separator_2, separator_3)
    ids, feature_vectors = domains_to_features(a_domains)

    results = {}
    if feature_vectors:
        classifier = load(MODEL)
        probabilities = classifier.predict_proba(feature_vectors)
        amino_acid_classes = classifier.classes_

        for i, seq_id in enumerate(ids):
            probability_list = probabilities[i]
            probs_and_aa = get_top_x_aa(amino_acid_classes, probability_list, nr_results)
            results[seq_id] = probs_and_aa

    else:
        print("No A domains found.")

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

    clear_temp()

    return results


def run_parasect(fasta_file, out_dir, job_name, extraction_method='hmm', save_signatures=False,
                 save_extended_signatures=False, save_domains=False, nr_results=3, separator_1=SEPARATOR_1,
                 separator_2=SEPARATOR_2, separator_3=SEPARATOR_3):

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    a_domains = get_domains(fasta_file, extraction_method, job_name, separator_1, separator_2, separator_3)
    ids, feature_vectors = domains_to_features(a_domains)

    results = {}
    if feature_vectors:
        classifier = load(MODEL)
        probabilities = classifier.predict_proba(feature_vectors)
        amino_acid_classes = classifier.classes_

        for i, seq_id in enumerate(ids):
            probability_list = probabilities[i]
            probs_and_aa = get_top_x_aa(amino_acid_classes, probability_list, nr_results)
            results[seq_id] = probs_and_aa

    else:
        print("No A domains found.")

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

    clear_temp()

    return results
