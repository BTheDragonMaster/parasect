from Bio import SeqIO
from collections import OrderedDict

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
        ID = None
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

        if ID is not None:
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

def proteins_from_genbank(gbk_in: str, fasta_out: str) -> None:
    """
    Convert a GBK file to a fasta file containing protein sequences

    Input:
    gbk_in: str, path to .gbk/.gb file
    fasta_out: str, path to output file

    """
    counter = 0
    id_to_seq = {}
    for record in SeqIO.parse(gbk_in, "genbank"):
        for feature in record.features:
            if feature.type == 'CDS':
                if 'translation' in feature.qualifiers:
                    sequence = feature.qualifiers['translation']

                    if 'protein_id' in feature.qualifiers:
                        seq_id = feature.qualifiers['protein_id'][0]
                    elif 'gene_id' in feature.qualifiers:
                        seq_id = feature.qualifiers['gene_id'][0]
                    elif 'locus_tag' in feature.qualifiers:
                        seq_id = feature.qualifiers['locus_tag'][0]
                    else:
                        counter += 1
                        seq_id = f"gene_{counter}"
                    id_to_seq[seq_id] = sequence

    write_fasta(id_to_seq, fasta_out)


class Tabular:
    def __init__(self, tabular_path, id_columns, separator='\t'):
        self.separator = separator
        self.data = OrderedDict()
        with open(tabular_path, 'r') as tabular_file:
            self.header = tabular_file.readline()
            header = self.header.split(self.separator)
            self.categories = []
            for i, category in enumerate(header):
                self.categories.append(category.strip())
            for line in tabular_file:

                values = line.split(self.separator)
                row_id = []
                for id_column in id_columns:
                    row_id.append(values[id_column])

                row_id = tuple(row_id)

                if row_id in self.data:
                    print(f"WARNING! Duplicate row ID found: {row_id}")
                self.data[row_id] = OrderedDict()

                for j, value in enumerate(values):
                    category = self.categories[j]
                    self.data[row_id][category] = value.strip()

    def get_column(self, category):
        if category not in self.categories:
            raise KeyError(f"Cannot find category {category} in data.")

        column = []

        for data_id in self.data:
            column.append(self.get_value(data_id, category))

        return column

    def get_row(self, data_id):
        if data_id not in self.data:
            raise KeyError(f"Cannot find data ID {data_id} in data.")

        row = []

        for category in self.categories:
            row.append(self.get_value(data_id, category))

        return row

    def get_value(self, data_id, category):
        try:
            return self.data[data_id][category]
        except KeyError:
            if data_id in self.data:
                print(f"Cannot find category {category} in data.")
            else:
                print(f"Cannot find id {data_id} in data.")
            raise KeyError

    def write_table(self, out_file):
        with open(out_file, 'w') as out:
            out.write(self.header)
            for seq_id in self.data:
                for i, category in enumerate(self.categories):
                    if i == len(self.categories) - 1:
                        out.write(f"{self.data[seq_id][category]}\n")
                    else:
                        out.write(f"{self.data[seq_id][category]}{self.separator}")


def parse_substrate_list(in_file):
    substrates = []
    with open(in_file, 'r') as substrate_file:
        for line in substrate_file:
            line = line.strip()
            if line:
                substrates.append(line)

    return substrates


def parse_amino_acid_properties(properties_file, return_categories=False):
    aa_to_vector = {}
    properties = Tabular(properties_file, [0])
    for data_id in properties.data:
        amino_acid = properties.get_value(data_id, 'AA')
        properties_string = properties.get_row(data_id)[1:]
        properties_float = list(map(float, properties_string))
        assert len(properties_float) == 15
        aa_to_vector[amino_acid] = properties_float

    if return_categories:
        return aa_to_vector, properties.categories[1:]

    return aa_to_vector


def parse_morgan_fingerprint(bitvector_file, return_categories=False):
    name_to_vector = OrderedDict()
    bitvector_data = Tabular(bitvector_file, [0])
    for name in bitvector_data.data:
        row = bitvector_data.get_row(name)
        name_to_vector[row[0]] = list(map(int, row[2:]))
    if return_categories:
        return name_to_vector, bitvector_data.categories[2:]
    return name_to_vector
