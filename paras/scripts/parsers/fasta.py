def read_fasta(fasta_dir):
    """
    Return dict of {seq_id: sequence, ->} from .fasta file

    Input:
    fasta_dir: str, directory of .fasta file

    Output:
    fasta_dict: dict of {seq_id: sequence, ->}, with seq_id and sequence str
    """
    with open(fasta_dir, 'r') as fasta_file:
        fasta_dict = {}
        sequence = []
        for line in fasta_file:
            line = line.strip()

            if line.startswith(">"):
                if sequence:
                    fasta_dict[ID] = ''.join(sequence)

                    ID = line[1:]
                    sequence = []
                else:
                    ID = line[1:]

            else:
                sequence.append(line)

        fasta_dict[ID] = ''.join(sequence)
        fasta_file.close()
    return fasta_dict


def write_fasta(id_to_seq, out_file):
    """Write .fasta file from dictionary containing sequence ids and sequences

    Input:
    id_to_seq: dict of {id: seq, ->}, with id and seq str
    out_file: str, file location of output .fasta file
    """
    with open(out_file, 'w') as out:
        for id, seq in id_to_seq.items():
            out.write(f'>{id}\n{seq}\n')
