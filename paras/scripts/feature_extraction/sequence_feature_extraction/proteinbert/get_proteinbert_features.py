from sys import argv
import os

from proteinbert import load_pretrained_model
from proteinbert.conv_and_global_attention_model import get_model_with_hidden_layers_as_outputs
from paras.scripts.parsers.fasta import read_fasta


seq_len = 1000
batch_size = 200


def get_feature_vectors(fasta_file, out_dir):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    seqid_to_fasta = read_fasta(fasta_file)
    seq_ids = []
    seqs = []
    for seq_id, seq in seqid_to_fasta.items():
        seq_ids.append(seq_id)
        seqs.append(seq)

    pretrained_model_generator, input_encoder = load_pretrained_model()
    model = get_model_with_hidden_layers_as_outputs(pretrained_model_generator.create_model(seq_len))
    encoded_x = input_encoder.encode_X(seqs, seq_len)
    local_representations, global_representations = model.predict(encoded_x, batch_size=batch_size)
    return seq_ids, local_representations, global_representations
