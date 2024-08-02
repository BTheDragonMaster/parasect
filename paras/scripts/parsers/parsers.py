from typing import List, Optional, Dict
from paras.scripts.parsers.fasta import read_fasta
from paras.scripts.parsers.datapoint import DataPoint
from paras.scripts.parsers.tabular import Tabular
from paras.scripts.parsers.voxel import Voxel, VoxelImportanceVector
from paras.scripts.math.shapes import Vector3D, Cube
from collections import OrderedDict
from dataclasses import dataclass, field


def parse_taxonomy(taxonomy_file):
    protein_to_taxonomy: dict[str, list[str]] = {}
    with open(taxonomy_file, 'r') as tax_info:
        for line in tax_info:
            tax_data = line.strip().split('\t')
            protein = tax_data[0]
            taxonomy = tax_data[1:]
            protein_to_taxonomy[protein] = taxonomy
    return protein_to_taxonomy


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


@dataclass
class SubstrateResult:
    rank: int
    substrate: str
    confidence: float


@dataclass
class ParasResult:
    protein_id: str
    domain_nr: int
    coordinates: tuple[int, int]
    signature: Optional[str] = None
    substrates: List[SubstrateResult] = field(default_factory=list)


def parse_paras_domain_name(domain_name: str, separator: str = '|') -> tuple[str, int, int, int]:
    """
    Return the domain name, number, and coordinates from a paras results domain identifier

    Parameters
    ----------
    domain_name: str, paras results domain identifier
    separator: str, separator used to make paras results

    Returns
    -------
    protein_id: str, identifier of the protein
    domain_nr: int, domain index
    start: int, start coordinate of domain within protein
    end: int, end coordinate of domain within protein

    """
    domain_info = domain_name.split(separator)
    protein_id: str = separator.join(domain_info[:-2])
    domain_nr: int = int(domain_info[-2].split('_')[1])
    coordinates: list[str] = domain_info[-1].split('-')
    start: int = int(coordinates[0])
    end: int = int(coordinates[1])

    return protein_id, domain_nr, start, end


def parse_results(results_file, signature_file: Optional[str] = None):
    results = Tabular(results_file, [0])
    domain_to_result = {}
    id_to_signature = {}
    if signature_file:
        id_to_signature = read_fasta(signature_file)
    for datapoint in results.data:
        domain_name = results.get_value(datapoint, "sequence_id")

        protein_id, domain_nr, start, end = parse_paras_domain_name(domain_name)
        coordinates = (start, end)

        if signature_file:
            signature = id_to_signature[domain_name]
            domain = ParasResult(protein_id, domain_nr, coordinates, signature)
        else:
            domain = ParasResult(protein_id, domain_nr, coordinates)
        for i, category in enumerate(results.categories):
            if 'confidence_score' in category:
                rank = int(category.split('_')[-1])
                substrate = results.get_value(datapoint, results.categories[i - 1])
                confidence = float(results.get_value(datapoint, category))
                substrate_result = SubstrateResult(rank, substrate, confidence)
                domain.substrates.append(substrate_result)
        domain.substrates.sort(key = lambda x: x.rank)
        domain_to_result[domain_name] = domain

    return domain_to_result


def yield_results(results_file):
    results = Tabular(results_file, [0])
    for datapoint in results.data:
        domain_name = results.get_value(datapoint, "sequence_id")
        protein_id, domain_nr, start, end = parse_paras_domain_name(domain_name)
        coordinates = (start, end)

        domain = ParasResult(protein_id, domain_nr, coordinates)
        for i, category in enumerate(results.categories):
            if 'confidence_score' in category:
                rank = int(category.split('_')[1])
                substrate = results.get_value(datapoint, results.categories[i - 1])
                confidence = float(results.get_value(datapoint, category))
                substrate_result = SubstrateResult(rank, substrate, confidence)
                domain.substrates.append(substrate_result)
        yield domain


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
        for i, amino_acid in enumerate(amino_acids[:]):
            amino_acids[i] = amino_acid.strip()
        for line in cm_data:
            row_repr = line.split('\t')[1:]
            row = []
            for element in row_repr:
                row.append(int(element))

            matrix.append(row)

        return matrix, amino_acids


@dataclass
class Prediction:
    predictor: str
    predicted_substrates: List[str]
    true_substrates: List[str]
    correct: bool
    no_call: bool


@dataclass
class SubstrateMetrics:
    substrate: str
    tp: int = 0
    fp: int = 0
    fn: int = 0

    precision: float = 0
    recall: float = 0
    f1: float = 0

    def set_metrics(self):
        if self.tp or self. fp:
            self.precision = float(self.tp) / (self.tp + self.fp)
        else:
            self.precision = 0.0
        if self.tp or self.fn:
            self.recall = float(self.tp) / (self.tp + self.fn)
        else:
            self.recall = 0.0

        if self.precision > 0.0 or self.recall > 0.0:
            self.f1 = 2 * self.precision * self.recall / (self.precision + self.recall)
        else:
            self.f1 = 0.0


@dataclass
class ToolPerformance:
    predictor: str
    predictions: List[Prediction] = field(default_factory=list)
    substrate_metrics: Dict[str, SubstrateMetrics] = field(default_factory=dict)

    def set_per_substrate_metrics(self, included_substrates):
        for substrate in included_substrates:
            substrate_metrics = SubstrateMetrics(substrate)

            for prediction in self.predictions:
                if prediction.correct and substrate in prediction.true_substrates:
                    substrate_metrics.tp += 1
                elif substrate in prediction.true_substrates:
                    substrate_metrics.fn += 1
                elif substrate in prediction.predicted_substrates and not prediction.correct:
                    substrate_metrics.fp += 1

            substrate_metrics.set_metrics()
            self.substrate_metrics[substrate] = substrate_metrics


def parse_summary_file(summary_file):
    summary = Tabular(summary_file, [0])
    predictor_to_performance = {}

    for datapoint in summary.data:
        true_substrates = summary.get_value(datapoint, "substrate").split('|')

        for predictor in summary.categories[2:]:

            if predictor not in predictor_to_performance:
                predictor_to_performance[predictor] = ToolPerformance(predictor)

            substrate_predictions = summary.get_value(datapoint, predictor).split('|')

            is_correct = False
            no_call = False
            for prediction in substrate_predictions:

                if prediction in true_substrates:
                    is_correct = True
                elif prediction in ['no_call', 'N/A', "no_confident_result", "no_force_needed"]:
                    no_call = True

            prediction = Prediction(predictor, substrate_predictions, true_substrates, is_correct, no_call)
            predictor_to_performance[predictor].predictions.append(prediction)

    return predictor_to_performance


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


def parse_sequence_feature_importances(importance_file, has_header=False):
    feature_to_importance = {}
    with open(importance_file, 'r') as importances:
        if has_header:
            importances.readline()
        for line in importances:
            line = line.strip()
            feature, importance = line.split('\t')
            feature_to_importance[feature] = float(importance)

    return feature_to_importance


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


def parse_proteinbert_features(proteinbert_file, return_categories=False):
    id_to_proteinbert = OrderedDict()
    feature_nr = None
    with open(proteinbert_file, 'r') as proteinbert:
        for line in proteinbert:
            line_data = line.strip().split('\t')
            seq_id = line_data[0]
            features = list(map(float, line_data[1:]))
            feature_nr = len(features)
            id_to_proteinbert[seq_id] = features

    if not return_categories:
        return id_to_proteinbert
    else:
        categories = []
        for i in range(feature_nr):
            categories.append(f"ProteinBERT_{i + 1}")
        return id_to_proteinbert, categories



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
