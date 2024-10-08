import os

from paras.parsers import proteins_from_genbank
from paras.feature_extraction import domains_from_fasta
from paras.feature_extraction import rename_sequences, reverse_renaming


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


