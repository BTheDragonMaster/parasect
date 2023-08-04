from paras.scripts.parsers.tabular import Tabular
from paras.scripts.parsers.parsers import yield_identities, parse_specificities, parse_domain_list
from paras.scripts.parsers.fasta import read_fasta
from sys import argv
import networkx as nx
import math
from pprint import pprint


def deduplicate_sequences(sequence_file, out_file):
    domains = []
    sequence_data = Tabular(sequence_file, [0])
    for domain_id in sequence_data.data:
        domain = Domain(domain_id)
        domain.add_sequence('mibig', sequence_data.get_value(domain_id, 'mibig_seq'))
        domain.add_specificity('mibig', sequence_data.get_value(domain_id, 'mibig_specificity'))
        domain.add_sequence('paras', sequence_data.get_value(domain_id, 'paras_seq'))
        domain.add_specificity('paras', sequence_data.get_value(domain_id, 'paras_specificity'))
        domain.add_sequence('sandpuma', sequence_data.get_value(domain_id, 'sandpuma_seq'))
        domain.add_specificity('sandpuma', sequence_data.get_value(domain_id, 'sandpuma_specificity'))
        domain.add_sequence('nrpspredictor', sequence_data.get_value(domain_id, 'nrpspredictor_seq'))
        domain.add_specificity('nrpspredictor', sequence_data.get_value(domain_id, 'nrpspredictor_specificity'))
        domain.set_consensus_sequence()
        domains.append(domain)

    G = nx.DiGraph()
    for domain in domains:
        G.add_node(domain)

    for domain_1 in domains:
        for domain_2 in domains:
            if domain_1.seq_id != domain_2.seq_id:
                if (domain_1.sequence[30:] in domain_2.sequence and domain_2.sequence[:100] in domain_1.sequence) or \
                        (domain_1.sequence[:100] in domain_2.sequence and domain_2.sequence[30:] in domain_1.sequence) or \
                        (domain_1.sequence in domain_2.sequence) or (domain_2.sequence in domain_1.sequence) or \
                        (domain_1.sequence == domain_2.sequence):
                    G.add_edge(domain_1, domain_2)
                    G.add_edge(domain_2, domain_1)

    connected_components = list(nx.strongly_connected_components(G))
    domain_groups = []
    for component in connected_components:
        domain_groups.append(DomainGroup(component))

    with open(out_file, 'w') as out:
        out.write(
            "domain_id\tsequence\tmibig_specificity\tparas_specificity\tsandpuma_specificity\tnrpspredictor_specificity\n")

        for domain in domain_groups:
            out.write(f"{domain.seq_id}\t{domain.sequence}\t{'|'.join(domain.mibig_spec)}\t{'|'.join(domain.paras_spec)}\t{'|'.join(domain.sandpuma_spec)}\t{'|'.join(domain.nrpspredictor_spec)}\n")


def has_specificity_conflict(domain_group):
    assessed_pairs = set()
    has_conflict = False
    for domain_1 in domain_group:
        for domain_2 in domain_group:
            pair = tuple(sorted([domain_1, domain_2], key=lambda x: x.seq_id))
            if domain_1 != domain_2 and pair not in assessed_pairs:
                if domain_1.mibig_spec and domain_2.mibig_spec and not domain_1.mibig_spec.intersection(domain_2.mibig_spec):
                    print(f"Conflict between mibig specificity of {domain_1} and {domain_2}:")
                    print(domain_1.mibig_spec)
                    print(domain_2.mibig_spec)
                    has_conflict = True
                if domain_1.paras_spec and domain_2.paras_spec and not domain_1.paras_spec.intersection(domain_2.paras_spec):
                    print(f"Conflict between paras specificity of {domain_1} and {domain_2}:")
                    print(domain_1.paras_spec)
                    print(domain_2.paras_spec)
                    has_conflict = True
                if domain_1.sandpuma_spec and domain_2.sandpuma_spec and not domain_1.sandpuma_spec.intersection(domain_2.sandpuma_spec):
                    print(f"Conflict between sandpuma specificity of {domain_1} and {domain_2}:")
                    print(domain_1.sandpuma_spec)
                    print(domain_2.sandpuma_spec)
                    has_conflict = True
                if domain_1.nrpspredictor_spec and domain_2.nrpspredictor_spec and not domain_1.nrpspredictor_spec.intersection(domain_2.nrpspredictor_spec):
                    print(f"Conflict between nrpspredictor specificity of {domain_1} and {domain_2}:")
                    print(domain_1.nrpspredictor_spec)
                    print(domain_2.nrpspredictor_spec)
                    has_conflict = True
                assessed_pairs.add(pair)
    return has_conflict


class DomainGroup:
    def __init__(self, domains):
        self.seq_id = '|'.join([domain.seq_id for domain in domains])

        self.domains = domains
        self.sequence = ''
        self.set_sequence()
        self.specificity = set()
        self.mibig_spec = set()
        self.paras_spec = set()
        self.sandpuma_spec = set()
        self.nrpspredictor_spec = set()
        self.set_specificities()


    def set_sequence(self):
        seq_length = 0
        for domain in self.domains:
            if len(domain.sequence) > seq_length:
                seq_length = len(domain.sequence)
                self.sequence = domain.sequence

    def set_specificities(self):
        assert not has_specificity_conflict(self.domains)

        for domain in self.domains:
            self.mibig_spec = self.mibig_spec.union(domain.mibig_spec)
            self.paras_spec = self.paras_spec.union(domain.paras_spec)
            self.sandpuma_spec = self.sandpuma_spec.union(domain.sandpuma_spec)
            self.nrpspredictor_spec = self.nrpspredictor_spec.union(domain.nrpspredictor_spec)


class Domain:

    def __init__(self, seq_id):
        self.seq_id = seq_id
        self.mibig_seq = ''
        self.paras_seq = ''
        self.sandpuma_seq = ''
        self.nrpspredictor_seq = ''
        self.sequence = ''
        self.specificity = set()

        self.mibig_spec = set()
        self.paras_spec = set()
        self.sandpuma_spec = set()
        self.nrpspredictor_spec = set()

    def __hash__(self):
        return hash(self.seq_id)

    def __eq__(self, other):
        return type(self) == type(other) and self.seq_id == other.seq_id

    def __repr__(self):
        return self.seq_id

    def add_sequence(self, label, sequence):
        if label == 'mibig':
            self.mibig_seq = sequence
        elif label == 'paras':
            self.paras_seq = sequence
        elif label == 'sandpuma':
            self.sandpuma_seq = sequence
        elif label == 'nrpspredictor':
            self.nrpspredictor_seq = sequence

    def add_specificity(self, label, specificity):
        if label == 'mibig':
            if specificity:
                self.mibig_spec = set(specificity.split('|'))
        elif label == 'paras':
            if specificity:
                self.paras_spec = set(specificity.split('|'))
        elif label == 'sandpuma':
            if specificity:
                self.sandpuma_spec = set(specificity.split('|'))
        elif label == 'nrpspredictor':
            if specificity:
                self.nrpspredictor_spec = set(specificity.split('|'))

    def has_sequence_conflicts(self):
        sequences = [self.mibig_seq, self.paras_seq, self.sandpuma_seq, self.nrpspredictor_seq]
        evaluated_pairs = set()
        problems_found = set()
        for i, seq_1 in enumerate(sequences):
            for j, seq_2 in enumerate(sequences):
                if i != j and seq_1 and seq_2 and tuple(sorted([i, j])) not in evaluated_pairs:
                    if (seq_1[30:] in seq_2 and seq_2[:100] in seq_1) or (seq_1[:100] in seq_2 and seq_2[30:] in seq_1) or (seq_1 in seq_2) or (seq_2 in seq_1) or (seq_1 == seq_2):
                        pass
                    else:
                        problems_found.add(tuple(sorted([i, j])))
                    evaluated_pairs.add(tuple(sorted([i, j])))

        labels = ["mibig", "paras", "sandpuma", "nrpspredictor"]

        if problems_found:
            print(f"Problems found with {self.seq_id}:")
            for problem_found in problems_found:
                print(f"Clash detected between {labels[problem_found[0]]} and {labels[problem_found[1]]}")
            return True
        return False

    def set_sequence(self, sequence):
        self.sequence = sequence

    def set_consensus_sequence(self):

        if self.has_sequence_conflicts():
            if self.mibig_seq:
                self.sequence = self.mibig_seq
            elif self.paras_seq:
                self.sequence = self.paras_seq
            elif self.sandpuma_seq:
                self.sequence = self.sandpuma_seq
            elif self.nrpspredictor_seq:
                self.sequence = self.nrpspredictor_seq
            else:
                raise ValueError(f"No sequence recorded for domain {self.seq_id}")
        else:
            seq_length = 0
            for sequence in [self.mibig_seq, self.paras_seq, self.sandpuma_seq, self.nrpspredictor_seq]:
                if len(sequence) > seq_length:
                    self.sequence = sequence
                    seq_length = len(sequence)

    def has_specificity_conflicts(self, exact=False):
        has_conflicts = False
        spec_groups = [self.mibig_spec, self.paras_spec, self.sandpuma_spec, self.nrpspredictor_spec]

        for i, spec_1 in enumerate(spec_groups):
            for j, spec_2 in enumerate(spec_groups):
                if j > i:
                    if exact:
                        if spec_1 and spec_2 and spec_1 != spec_2:
                            has_conflicts = True
                    else:
                        if spec_1 and spec_2 and not spec_1.intersection(spec_2):
                            has_conflicts = True
        return has_conflicts

    def set_consensus_specificity(self):
        assert not self.has_specificity_conflicts()
        spec_groups = [self.mibig_spec, self.paras_spec, self.sandpuma_spec, self.nrpspredictor_spec]
        for spec_group in spec_groups:
            if spec_group:
                self.specificity = spec_group
                break

    def print_specificities(self):
        print(self.seq_id)
        print("MIBiG:", '|'.join(self.mibig_spec))
        print("PARAS:", '|'.join(self.paras_spec))
        print("SANDPUMA:", '|'.join(self.sandpuma_spec))
        print("NRPSPredictor:", '|'.join(self.nrpspredictor_spec))
        print('\n')


def find_duplicates_from_identity(identities_file):
    domain_to_identical = {}
    for domain_1, domain_2, identity in yield_identities(identities_file):
        if math.isclose(100.0, identity):
            if domain_1 not in domain_to_identical:
                domain_to_identical[domain_1] = []
            domain_to_identical[domain_1].append(domain_2)
    pprint(domain_to_identical)


def find_clashes(duplicates_file, specificities_file, fasta_file, domain_list, out_file):
    domain_to_spec = parse_specificities(specificities_file)
    domain_to_seq = read_fasta(fasta_file)
    domains_to_remove = []
    for domain_1, domain_2, identity in yield_identities(duplicates_file):
        if math.isclose(100.0, identity):
            spec_1 = domain_to_spec[domain_1]
            spec_2 = domain_to_spec[domain_2]
            seq_1 = domain_to_seq[domain_1]
            seq_2 = domain_to_seq[domain_2]
            if not set(spec_1).intersection(set(spec_2)):
                print("Clash detected:", domain_1, spec_1, domain_2, spec_2)
            else:
                domain_to_keep = domain_1
                domain_to_remove = domain_2
                if len(seq_2) > len(seq_1):
                    domain_to_keep = domain_2
                    domain_to_remove = domain_1

                spec_to_keep = spec_1
                if len(spec_2) > len(spec_1):
                    spec_to_keep = spec_2

                if domain_to_remove in domains_to_remove:
                    print(domain_to_remove, spec_1, spec_2)

                domains_to_remove.append(domain_to_remove)

                print(domain_to_keep, spec_to_keep)

        print(set(domains_to_remove))

        with open(out_file, 'w') as out:
            for domain in parse_domain_list(domain_list):
                if domain not in domains_to_remove:
                    out.write(f"{domain}\n")

def filter_duplicates_from_list(identities_file, fasta_file, domain_list, out_file):
    domain_to_seq = read_fasta(fasta_file)
    domains_to_remove = []
    for domain_1, domain_2, identity in yield_identities(identities_file):
        if math.isclose(100.0, identity):

            seq_1 = domain_to_seq[domain_1]
            seq_2 = domain_to_seq[domain_2]

            domain_to_remove = domain_2
            if len(seq_2) > len(seq_1):

                domain_to_remove = domain_1

            domains_to_remove.append(domain_to_remove)

    print(set(domains_to_remove))

    with open(out_file, 'w') as out:
        for domain in parse_domain_list(domain_list):
            if domain not in domains_to_remove:
                out.write(f"{domain}\n")


if __name__ == "__main__":
    # deduplicate_sequences(argv[1], argv[2])
    # find_duplicates_from_identity(argv[1])
    # find_clashes(argv[1], argv[2], argv[3], argv[4], argv[5])
    filter_duplicates_from_list(argv[1], argv[2], argv[3], argv[4])

