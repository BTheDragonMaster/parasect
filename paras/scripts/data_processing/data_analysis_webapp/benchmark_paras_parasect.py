from sqlalchemy import create_engine
from sqlalchemy.orm import Session
from sys import argv
from typing import Optional, List
from dataclasses import dataclass, field

from parasect.database.query_database import get_domains_from_synonym
from parasect.core.parsing import parse_fasta_file, parse_list
from parasect.core.tabular import Tabular


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
    results = Tabular(results_file)
    domain_to_result = {}
    id_to_signature = {}
    if signature_file:
        id_to_signature = parse_fasta_file(signature_file)
    for datapoint in results.rows:
        domain_name = results.get_row_value(datapoint, "sequence_id")

        protein_id, domain_nr, start, end = parse_paras_domain_name(domain_name)
        coordinates = (start, end)

        if signature_file:
            signature = id_to_signature[domain_name]
            domain = ParasResult(protein_id, domain_nr, coordinates, signature)
        else:
            domain = ParasResult(protein_id, domain_nr, coordinates)
        for i, category in enumerate(results.column_names):
            if 'confidence_score' in category:
                rank = int(category.split('_')[-1])
                substrate = results.get_row_value(datapoint, results.column_names[i - 1])
                confidence = float(results.get_row_value(datapoint, category))
                substrate_result = SubstrateResult(rank, substrate, confidence)
                domain.substrates.append(substrate_result)
        domain.substrates.sort(key = lambda x: x.rank)
        domain_to_result[domain_name] = domain

    return domain_to_result


if __name__ == "__main__":
    engine = create_engine(f"sqlite:///{argv[1]}")
    paras_results = parse_results(argv[2])
    train_domains = parse_list(argv[3])
    correct = 0
    incorrect = 0
    in_training_data = 0
    with Session(engine) as session:
        for seq_id, result in paras_results.items():
            synonyms = seq_id.split('|')
            in_dataset = False
            for synonym in synonyms:
                for train_domain in train_domains:
                    if synonym in train_domain.split('|'):
                        in_dataset = True
            if in_dataset:
                in_training_data += 1
            else:
                synonym = synonyms[0]
                domain = get_domains_from_synonym(session, synonym)[0]
                true_substrates = [s.name for s in domain.substrates]
                predicted_substrate = result.substrates[0].substrate
                if predicted_substrate in true_substrates:
                    correct += 1
                else:
                    incorrect += 1

    print(in_training_data)
    print(correct, incorrect)
    print(correct / (incorrect + correct))


