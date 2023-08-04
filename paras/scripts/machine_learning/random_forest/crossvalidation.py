from sklearn.ensemble import RandomForestClassifier
from argparse import ArgumentParser
from paras.scripts.feature_extraction.sequence_feature_extraction.seq_to_features import get_sequence_features, PROPERTIES_FILE


from paras.scripts.parsers.parsers import parse_moment_vectors, parse_bitvectors, parse_specificities, \
    parse_domain_list, parse_pca_data
from paras.scripts.parsers.fasta import read_fasta
from paras.scripts.parsers.tabular import Tabular


def parse_arguments():
    parser = ArgumentParser()
    parser.add_argument('-i', "--input_folder", type=str, required=True, help="Path to crossvalidation directory.")
    parser.add_argument('-o', "--output_folder", type=str, required=True, help="Path to output directory.")
    featurisation = parser.add_mutually_exclusive_group(required=True)
    featurisation.add_argument('-pca', action="store_true", help="Featurise domains using x principal components of voxel data. Default x is 10; change by changing parameter -p.")
    featurisation.add_argument('-seq', action="store_true",
                               help="Featurise domains using 34 amino acid signature extracted using a structure-based HMM.")
    featurisation.add_argument('-mom', action="store_true",
                               help="Featurise domains using moment invariants for atom coordinates per atom type.")
    featurisation.add_argument('-pca_seq', action="store_true",
                               help="combines 'pca' and 'seq'.")
    featurisation.add_argument('-mom_seq', action="store_true", help="combines 'mom' and 'seq'.")

    model = parser.add_argument_group(required=True)
    model.add_argument("-ml", "--multi-label", action="store_true",
                       help="Build one model with multiple-label output. Only domain features are considered.")
    model.add_argument("-sl", "--single-label", action="store_true",
                       help="Build one model with single-label output. Only domain features are considered. Warning: only the first listed substrate is used.")
    model.add_argument("-bi", "--binary", action="store_true",
                       help="Build one model with binary output that featurises both domain and compound in tandem.")
    model.add_argument("-mm", "--multi-model", action="store_true",
                       help="Build one model per substrate with binary output.")
    parser.add_argument('-cf', "--compound_features", type=str, help="Directory to file containing bitvector data. Required if model type is 'binary'.")
    parser.add_argument('-df')

    parser.add_argument('-n', type=int, default=100, help="Nr of trees in RF")
    parser.add_argument('-t', type=int, default=None, help="Nr of threads. -1 uses all available cores")
    parser.add_argument('-sampling', type=str, default=None, help="balanced, balanced_subsample, over_sample, under_sample or None")

    args = parser.parse_args()
    return args


class RandomForestDataset:
    def __init__(self, domains, specificities_file, model_type, featurisation_method, feature_file):
        self.domains = domains
        self.model_type = model_type
        self.featurisation_method = featurisation_method
        self.substrates, self.all_substrates = self.get_substrates(specificities_file)
        self.sequence_feature_names = Tabular(PROPERTIES_FILE, [0]).categories[1:]


    def get_domain_feature_vectors(self, feature_file):
        feature_vectors = []
        if self.featurisation_method == 'pca':
            domain_to_vector = parse_pca_data(feature_file)
        elif self.featurisation_method == 'seq':
            domain_to_vector = self.get_domain_feature_vectors(feature_file)
        elif self.featurisation_method == 'pca_seq':
            domain_to_vector = {}
            domain_to_pca = parse_pca_data(feature_file)


    def get_sequence_features(self, feature_file):
        domain_to_seq = read_fasta(feature_file)
        domain_to_vector = {}
        for domain, seq in domain_to_seq.items():
            sequence_vector = get_sequence_features(seq)
            domain_to_vector[domain] = sequence_vector
        return domain_to_vector

    def get_substrates(self, specificities_file):
        specificities = []
        domain_to_specificities = parse_specificities(specificities_file)
        all_specificities = []
        for specificity_list in domain_to_specificities.values():
            for specificity in specificity_list:
                if specificity not in all_specificities:
                    all_specificities.append(specificity)
        all_specificities.sort()

        for domain in self.domains:
            specificity = domain_to_specificities[domain]

            if self.model_type == 'single-label':
                specificities.append(specificity[0])
            elif self.model_type == 'multi-model':
                specificity_vector = []
                for spec in all_specificities:
                    if spec in specificity:
                        specificity_vector.append(1)
                    else:
                        specificity_vector.append(0)
            else:
                specificities.append(specificity)

        return specificities, all_specificities


    def get_domain_featurisation(self, domain):



    def get_compound_featurisation(self, bitvector_file):



class ParasectRandomForest:
    def __init__(self, featurisation_method, model_type, nr_trees, sampling):
        self.featurisation_method = featurisation_method
        self.model_type = model_type
        self.tree_nr = nr_trees
        self.sampling_method = sampling

    def train(self, train_domains, domain_featurisations=None, compound_featurisations=None):
