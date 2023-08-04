import os

from paras.scripts.parsers.fasta import read_fasta
from paras.scripts.parsers.datapoint import DataPoint
from paras.scripts.parsers.tabular import Tabular
from paras.scripts.parsers.voxel import Voxel, VoxelImportanceVector
from paras.scripts.math.shapes import Vector3D, Cube
from collections import OrderedDict


def parse_closest_identities(identities_file):
    domain_to_id = {}
    with open(identities_file, 'r') as identities:
        for line in identities:
            line = line.strip()
            if line:
                domain, identity = line.split('\t')
                identity = float(identity)
                domain_to_id[domain] = identity

    return domain_to_id

def yield_identities(identities_file):
    with open(identities_file, 'r') as identities:
        for line in identities:
            line = line.strip()
            if line:
                domain_1, domain_2, identity = line.split('\t')
                identity = float(identity)
                yield domain_1, domain_2, identity


def parse_mmseq_clusters(cluster_file):
    cluster_to_domain = {}
    with open(cluster_file, 'r') as cluster_data:
        for line in cluster_data:
            line = line.strip()
            if line:
                cluster, domain = line.split('\t')
                if cluster not in cluster_to_domain:
                    cluster_to_domain[cluster] = []
                cluster_to_domain[cluster].append(domain)

    return cluster_to_domain


def parse_test_results(test_results_file, get_confidence=False):
    test_results = Tabular(test_results_file, [0])

    domain_to_result = {}

    for datapoint in test_results.data:
        domain = test_results.get_value(datapoint, "domain")
        prediction = test_results.get_value(datapoint, "prediction")
        confidence = test_results.get_value(datapoint, "confidence")
        if get_confidence:
            domain_to_result[domain] = (prediction, float(confidence))
        else:
            domain_to_result[domain] = prediction

    return domain_to_result


def parse_cm_matrix(cm_file):

    with open(cm_file, 'r') as cm_data:
        matrix = []
        amino_acids = cm_data.readline().split('\t')[1:]
        for line in cm_data:
            row_repr = line.split('\t')[1:]
            row = []
            for element in row_repr:
                row.append(int(element))

            matrix.append(row)

        return matrix, amino_acids


def parse_substrate_list(in_file):
    substrates = []
    with open(in_file, 'r') as substrate_file:
        for line in substrate_file:
            line = line.strip()
            if line:
                substrates.append(line)

    return substrates

def parse_stach_codes(in_file):
    seq_data = Tabular(in_file, [0])
    id_to_stach = {}
    id_to_34 = {}
    for data_id in seq_data.data:
        seq_id = seq_data.get_value(data_id, 'domain_name')
        stach = seq_data.get_value(data_id, 'stachelhaus')
        seq_34 = seq_data.get_value(data_id, 'active_site')
        id_to_stach[seq_id] = stach
        id_to_34[seq_id] = seq_34

    return id_to_stach, id_to_34


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


def parse_feature_importances(importance_file, voxel_size=1.2):
    voxel_to_importances = OrderedDict()

    with open(importance_file, 'r') as importances:
        for line in importances:
            line = line.strip()
            voxel_info, importance = line.split('\t')
            coord_string, category = voxel_info.split('|')
            coords = list(map(float, coord_string.split('_')))
            cube = Cube(Vector3D(coords[0], coords[1], coords[2]), voxel_size)
            voxel = Voxel(cube, [], importances=VoxelImportanceVector())
            if voxel not in voxel_to_importances:
                voxel_to_importances[voxel] = VoxelImportanceVector()

            voxel_to_importances[voxel].update_vector(category, float(importance))

    for voxel, vector in voxel_to_importances.items():
        voxel.importances = vector

    return list(voxel_to_importances.keys())

def parse_bitscores(bitscore_file, domain_list):

    domain_to_domain_to_bitscore = {}
    with open(bitscore_file, 'r') as bitscore_info:
        for line in bitscore_info:
            line = line.strip()
            domain_1, domain_2, bitscore = line.split('\t')
            if domain_1 in domain_list and domain_2 in domain_list:
                bitscore = float(bitscore)
                if domain_1 not in domain_to_domain_to_bitscore:
                    domain_to_domain_to_bitscore[domain_1] = {}
                domain_to_domain_to_bitscore[domain_1][domain_2] = bitscore
                if domain_2 not in domain_to_domain_to_bitscore:
                    domain_to_domain_to_bitscore[domain_2] = {}
                domain_to_domain_to_bitscore[domain_2][domain_1] = bitscore
    domain_to_vector = {}
    for domain_1 in domain_list:
        vector = []
        for domain_2 in domain_list:
            if domain_1 == domain_2:
                vector.append(1000)
            else:
                vector.append(domain_to_domain_to_bitscore[domain_1][domain_2])
        domain_to_vector[domain_1] = vector

    return domain_to_vector


def parse_pocket_features(feature_file, return_categories=False):
    domain_to_pocket_vector = OrderedDict()
    pocket_data = Tabular(feature_file, [0])

    for name in pocket_data.data:
        row = pocket_data.get_row(name)
        domain_to_pocket_vector[row[0]] = list(map(int, row[1:]))
    if not return_categories:
        return domain_to_pocket_vector
    else:
        return domain_to_pocket_vector, pocket_data.categories[1:]


def yield_pocket_features(feature_file):
    with open(feature_file, 'r') as pocket_data:
        categories = pocket_data.readline()
        categories = categories.strip().split('\t')[1:]
        for line in pocket_data:
            domain_and_vector = line.strip().split('\t')
            domain = domain_and_vector[0]
            vector = list(map(int, domain_and_vector[1:]))
            yield domain, vector, categories


def parse_interactions(substrate_file):
    domain_to_interaction = OrderedDict()
    with open(substrate_file, 'r') as substrates:
        substrates.readline()
        for line in substrates:
            line = line.strip()
            if line:
                domain, interaction = line.split('\t')
                interaction = int(interaction)
                domain_to_interaction[domain] = interaction
    return domain_to_interaction


def parse_domain_list(domain_file):
    domain_ids = []
    with open(domain_file, 'r') as domains:
        for line in domains:
            line = line.strip()
            if line:
                domain_ids.append(line)
    return domain_ids


def parse_specificities(data_file):
    domain_to_spec = OrderedDict()
    parasect_data = Tabular(data_file, [0])
    for data_id in parasect_data.data:
        domain_id = parasect_data.get_value(data_id, "domain_id")
        domain_to_spec[domain_id] = parasect_data.get_value(data_id, "specificity").split('|')
    return domain_to_spec


def parse_pca_file(pca_file, return_categories=False):
    pca_data = Tabular(pca_file, [0])
    domain_to_vector = {}

    for domain in pca_data.data:
        domain_name = pca_data.get_value(domain, 'domain')
        vector = list(map(float, pca_data.get_row(domain)[2:]))
        domain_to_vector[domain_name] = vector
    if return_categories:
        return domain_to_vector, pca_data.categories[2:]
    return domain_to_vector


def parse_pca_data(pca_file, return_categories=False):
    pca_data = Tabular(pca_file, [0])
    domains = []
    pca_vectors = []

    for domain in pca_data.data:
        domains.append(domain)
        pca_vectors.append(list(map(float, pca_data.get_row(domain)[2:])))
    if return_categories:
        return pca_vectors, pca_data.categories[2:]
    return pca_vectors


def parse_metrics(metrics_file):
    metric_to_value = {}
    with open(metrics_file, 'r') as metrics:
        for line in metrics:
            line = line.strip()
            if line:
                metric_type, value = line.split('\t')
                if value == 'N/A':
                    pass
                elif metric_type in ["correct", "incorrect"]:
                    value = int(value)
                else:
                    value = float(value)

                metric_to_value[metric_type] = value

    return metric_to_value


def parse_moment_vectors(moment_file):
    name_to_vector = OrderedDict()
    moment_data = Tabular(moment_file, [0])
    for name in moment_data.data:
        row = moment_data.get_row(name)
        name_to_vector[row[0]] = list(map(float, row[1:]))
    return name_to_vector


def parse_morgan_fingerprint(bitvector_file, return_categories=False):
    name_to_vector = OrderedDict()
    bitvector_data = Tabular(bitvector_file, [0])
    for name in bitvector_data.data:
        row = bitvector_data.get_row(name)
        name_to_vector[row[0]] = list(map(int, row[2:]))
    if return_categories:
        return name_to_vector, bitvector_data.categories[2:]
    return name_to_vector


def parse_substrate_smiles(smiles_file):
    substrate_to_smiles = {}
    smiles_data = Tabular(smiles_file, [0])
    for name in smiles_data.data:
        substrate_id = smiles_data.get_value(name, "substrate")
        smiles_string = smiles_data.get_value(name, "smiles")
        substrate_to_smiles[substrate_id] = smiles_string

    return substrate_to_smiles


def parse_paras_specificities(specificities_file):
    domain_to_specificities = {}
    with open(specificities_file, 'r') as specificities:
        for line in specificities:
            line = line.strip()
            if line:
                domain, specificity = line.split('\t')
                domain_to_specificities[domain] = specificity.split('|')

    return domain_to_specificities


def parse_mibig_specificities(specificities_file):
    domain_to_data = {}
    mibig_data = Tabular(specificities_file, [2, 4])
    for seq_id in mibig_data.data:
        protein_id = mibig_data.data[seq_id]["protein_name"]
        a_domain_nr = mibig_data.data[seq_id]["domain_nr"]
        datapoint = DataPoint()
        datapoint.set_id(f"{protein_id}.A{a_domain_nr}")
        if mibig_data.data[seq_id]['sequence']:
            datapoint.set_sequence(mibig_data.data[seq_id]['sequence'])
        else:
            pass


def parse_sandpuma_dataset(sandpuma_fasta, abbreviations_file):
    header_to_sequence = read_fasta(sandpuma_fasta)

    abbreviations = Tabular(abbreviations_file, [0])

    id_to_specificity = {}

    for header, sequence in header_to_sequence.items():
        adom_id, specificity_abbr, protein_id, adom_nr, bgc_id, compound = header.split('\t')
        specificity_abbreviations = specificity_abbr.split('|')
        specificities = []
        for abbreviation in specificity_abbreviations:
            full_name = abbreviations.get_value(abbreviation, "Verbose name")
            specificities.append(full_name)

        seq_id = f"{protein_id}.{adom_nr}"
        id_to_specificity[seq_id] = specificities

    return id_to_specificity
