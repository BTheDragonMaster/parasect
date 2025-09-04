from argparse import ArgumentParser, Namespace
import os
from enum import Enum

from parasect.database.query_database import sequences_are_equivalent
from parasect.core.parsing import parse_parasect_data


def parse_arguments() -> Namespace:
    """Parse command line arguments

    :return: Arguments
    :rtype: Namespace
    """

    parser = ArgumentParser(description="Compare two training datasets")
    parser.add_argument('-t1', type=str, required=True,
                        help="Path to first training set in parasect dataset format")
    parser.add_argument('-t2', type=str, required=True,
                        help="Path to second training set in parasect dataset format")
    parser.add_argument('-o', type=str, required=True,
                        help="Path to output folder")
    parser.add_argument('-l1', type=str, default="t1",
                        help="Label for first training set")
    parser.add_argument('-l2', type=str, default="t2",
                        help="Label for second training set")

    arguments = parser.parse_args()

    return arguments


class SpecificityMatchType(Enum):
    MATCH = 1
    MISMATCHING_ORDER = 2
    PARTIAL_MISMATCH = 3
    FULL_MISMATCH = 4


def compare_specificities(specificities_1: list[str], specificities_2: list[str]) -> SpecificityMatchType:
    """

    :param specificities_1: list of specificities for a single domain
    :type specificities_1: list[str]
    :param specificities_2: list of specificities for a single domain
    :type specificities_2: list[str]
    :return: Type of match/mismatch
    :rtype: SpecificityMatchType
    """

    if specificities_1 != specificities_2:
        if set(specificities_1) == set(specificities_2):
            return SpecificityMatchType.MISMATCHING_ORDER
        elif set(specificities_1).intersection(set(specificities_2)):
            return SpecificityMatchType.PARTIAL_MISMATCH
        else:
            return SpecificityMatchType.FULL_MISMATCH
    else:
        return SpecificityMatchType.MATCH


def compare_training_sets(training_1: str, training_2: str,
                          label_1: str, label_2: str,
                          out_dir: str) -> None:
    """Compare sequences and specificities of two training sets

    :param training_1: Path to file containing first training set in parasect format
    :type training_1: str
    :param training_2: Path to file containing second training set in parasect format
    :type training_2: str
    :param label_1: Label of first training set
    :type label_1: str
    :param label_2: Label of second training set
    :type label_2: str
    :param out_dir: Path to output directory. Will contain mapping file, specificity mismatch file,
        and files containing domains unique to each dataset
    :type out_dir: str
    :return:
    :rtype:
    """

    id_to_seq_1, id_to_spec_1 = parse_parasect_data(training_1)
    id_to_seq_2, id_to_spec_2 = parse_parasect_data(training_2)

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    t1_exclusive_out = os.path.join(out_dir, f"{label_1}_exclusive.txt")
    t2_exclusive_out = os.path.join(out_dir, f"{label_2}_exclusive.txt")
    mismatches_out = os.path.join(out_dir, "mismatches.txt")
    id_mapping = os.path.join(out_dir, "id_mapping.txt")

    with open(t1_exclusive_out, 'w') as t1_out, open(t2_exclusive_out, 'w') as t2_out, \
            open(mismatches_out, 'w') as mis_out, open(id_mapping, 'w') as mapping_out:

        t2_ids = set(id_to_seq_2.keys())
        t1_to_t2: dict[str, list[str]] = {}
        mapping_out.write(f"{label_1}_id\t{label_2}_ids\n")
        mis_out.write(f"{label_1}_id\t{label_2}_id\tmismatch_type\tspecificities_{label_1}\tspecificities_{label_2}\n")
        for seq_id_1, seq_1 in id_to_seq_1.items():
            match_found = False
            for seq_id_2, seq_2 in id_to_seq_2.items():
                if sequences_are_equivalent(seq_1, seq_2):
                    match_found = True
                    if seq_id_1 not in t1_to_t2:
                        t1_to_t2[seq_id_1] = []
                    t1_to_t2[seq_id_1].append(seq_id_2)
                    if seq_id_2 in t2_ids:
                        t2_ids.remove(seq_id_2)

                    specs_1 = id_to_spec_1[seq_id_1]
                    specs_2 = id_to_spec_2[seq_id_2]

                    substrate_match = compare_specificities(specs_1, specs_2)
                    if substrate_match != SpecificityMatchType.MATCH:
                        mis_out.write(f"{seq_id_1}\t{seq_id_2}\t{substrate_match.name}\t{'|'.join(specs_1)}\t{'|'.join(specs_2)}\n")
            if not match_found:
                t1_out.write(f"{seq_id_1}\n")

        for seq_id in t2_ids:
            t2_out.write(f"{seq_id}\n")

        for seq_id_1, seq_ids_2 in t1_to_t2.items():
            mapping_out.write(f"{seq_id_1}\t{'|||'.join(seq_ids_2)}\n")


if __name__ == "__main__":
    args = parse_arguments()
    compare_training_sets(args.t1, args.t2, args.l1, args.l2, args.o)
