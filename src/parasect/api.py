import os
from collections import OrderedDict

from parasect.core.featurisation import domains_to_features, get_domains
from parasect.core.helpers import clear_temp_dir


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


def run_paras(
    selected_input: str,
    selected_input_type: str,
    temp_dir: str,
    first_separator: str,
    second_separator: str,
    third_separator: str,
    model,
    use_structure_guided_alignment: bool = False,
    num_predictions_to_report: int = 3,
    save_active_site_signatures: bool = False,
    save_extended_signatures: bool = False,
    save_adenylation_domain_sequences: bool = False,
):
    # write selected_input to file in temp folder
    file_name = "input.fasta" if selected_input_type == "fasta" else "input.gbk"
    input_file = os.path.join(temp_dir, file_name)
    with open(input_file, "w") as f:
        f.write(selected_input)

    # get domains
    a_domains = get_domains(
        input_file=input_file,
        extraction_method="profile" if use_structure_guided_alignment else "hmm",
        job_name="run",
        separator_1=first_separator,
        separator_2=second_separator,
        separator_3=third_separator,
        verbose=False,
        file_type=selected_input_type.lower(),
        temp_dir=temp_dir,
    )
    ids, feature_vectors = domains_to_features(a_domains, one_hot=False)

    # run model and retrieve class predictions
    results = OrderedDict()

    probabilities = model.predict_proba(feature_vectors)
    amino_acid_classes = model.classes_

    for i, seq_id in enumerate(ids):
        probability_list = probabilities[i]
        probs_and_aa = get_top_n_aa_paras(
            amino_acid_classes, probability_list, num_predictions_to_report
        )
        results[seq_id] = probs_and_aa

    model = None
    clear_temp_dir(temp_dir, keep=[".gitkeep"])

    # parse results
    domain_results = {}
    for domain in a_domains:
        domain_results[domain.domain_id] = {}
        if save_adenylation_domain_sequences:
            domain_results[domain.domain_id]["sequence"] = domain.sequence
        if save_active_site_signatures:
            domain_results[domain.domain_id]["signature"] = domain.signature
        if save_extended_signatures:
            domain_results[domain.domain_id]["extended_signature"] = domain.extended_signature

    for domain_id in results:
        preds = results[domain_id]
        domain_results[domain_id]["predictions"] = preds

    return [
        {"domain_id": domain_id, "data": domain_results[domain_id]} for domain_id in domain_results
    ]


def run_parasect(
    model,
    input_file,
    selected_input_type,
    first_separator,
    second_separator,
    third_separator,
    use_structure_guided_alignment,
    fingerprints,
    temp_dir,
    substrates,
    num_predictions_to_report=3,
    save_active_site_signatures=False,
    save_extended_signatures=False,
    save_adenylation_domain_sequences=False,
):
    a_domains = get_domains(
        input_file=input_file,
        extraction_method="profile" if use_structure_guided_alignment else "hmm",
        job_name="run",
        separator_1=first_separator,
        separator_2=second_separator,
        separator_3=third_separator,
        verbose=False,
        file_type=selected_input_type.lower(),
        temp_dir=temp_dir,
    )
    sequence_ids, sequence_feature_vectors = domains_to_features(a_domains, one_hot=False)

    if not sequence_feature_vectors:
        raise Exception("No feature vectors.")

    # Run model and retrieve class predictions.
    results = OrderedDict()

    batch_size = 1000
    counter = 0
    start = 0
    end = batch_size

    id_to_probabilities = {}

    batch_nr = 1
    while start < len(sequence_feature_vectors):

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

        probabilities = model.predict_proba(feature_vectors)
        interaction_labels = model.classes_

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
        results[seq_id] = get_top_n_aa_parasect(
            seq_id, id_to_probabilities, num_predictions_to_report
        )

    model = None
    clear_temp_dir(temp_dir, keep=[".gitkeep"])

    domain_results = {}
    for domain in a_domains:
        domain_results[domain.domain_id] = {}
        if save_adenylation_domain_sequences:
            domain_results[domain.domain_id]["sequence"] = domain.sequence
        if save_active_site_signatures:
            domain_results[domain.domain_id]["signature"] = domain.signature
        if save_extended_signatures:
            domain_results[domain.domain_id]["extended_signature"] = domain.extended_signature

    for domain_id in results:
        preds = results[domain_id]
        domain_results[domain_id]["predictions"] = preds

    return [
        {"domain_id": domain_id, "data": domain_results[domain_id]} for domain_id in domain_results
    ]
