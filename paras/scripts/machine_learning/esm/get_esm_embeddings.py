import os
from sys import argv

os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
os.environ["OMP_NUM_THREADS"] = "1"

import esm
import torch
from parasect.core.parsing import parse_fasta_file


def parse_positions_file(pos_file):
    domain_to_positions = {}
    with open(pos_file, 'r') as positions:
        positions.readline()
        for line in positions:
            line = line.strip()
            if line:
                position_data = line.split('\t')
                domain_name = position_data[0]
                signature_positions = position_data[1:]
                assert len(signature_positions) == 34
                if domain_name in domain_to_positions:
                    raise ValueError(f"Multiple entries for domain {domain_name}")
                domain_to_positions[domain_name] = []

                for position in signature_positions:
                    if position == 'None':
                        domain_to_positions[domain_name].append(None)
                    else:
                        domain_to_positions[domain_name].append(int(position))

    return domain_to_positions


if __name__ == "__main__":
    print("Loading ESM model..")
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model.eval()
    print("Preparing data..")
    batch_converter = alphabet.get_batch_converter()

    id_to_seq = parse_fasta_file(argv[1])
    new_id_to_seq = {}

    out_dir = argv[3]
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    for seq_id, seq in id_to_seq.items():
        safe_seq_id = seq_id.replace(r'/', '*')
        out_file = os.path.join(out_dir, f"{safe_seq_id}_esm.tsv")
        if not os.path.exists(out_file):
            new_id_to_seq[seq_id] = seq
        else:
            print(f"Already computed embeddings for {seq_id}. Skipping..")

    id_to_pos = parse_positions_file(argv[2])
    data = sorted(new_id_to_seq.items(), key=lambda x: len(x[1]))

    batch_labels, batch_strs, batch_tokens = batch_converter(data)

    print("Extracting per-residue tokens..")

    residue_embeddings_all = {}
    batch_size = 4

    for i in range(0, len(data), batch_size):
        batch = data[i:i + batch_size]
        batch_labels, batch_strs, batch_tokens = batch_converter(batch)

        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[33], return_contacts=False)
            token_representations = results["representations"][33]  # (B, L+2, 1280)

        for j, (seq_id, seq) in enumerate(batch):
            # Slice to get only actual residues (exclude [CLS] and [EOS])
            residue_representations = token_representations[j, 1:len(seq) + 1]
            positions = id_to_pos[seq_id]
            safe_seq_id = seq_id.replace(r'/', '*')
            out_file = os.path.join(out_dir, f"{safe_seq_id}_esm.tsv")
            with open(out_file, 'w') as out:
                out.write("residue")
                for k in range(1280):
                    out.write(f"\tesm_{k + 1}")
                out.write("\n")

                for k, position in enumerate(positions):
                    out.write(f"{k + 1}")
                    if position is None:
                        out.write("\t" + "\t".join(["0.0"] * 1280))
                    else:
                        out.write("\t" + "\t".join(map(str, residue_representations[position].tolist())))
                    out.write("\n")
            print(f"{seq_id}: {residue_representations.shape}")

