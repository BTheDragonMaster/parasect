import random
import os
from typing import Optional, Iterable
from argparse import ArgumentParser, Namespace
from enum import Enum

from sqlalchemy.orm import Session
from sqlalchemy import create_engine

from parasect.database.build_database import AdenylationDomain, Substrate
from parasect.core.taxonomy import Rank

from .domain_scope import DomainScope


def parse_arguments() -> Namespace:
    """Parse arguments from command line

    :return: Arguments
    :rtype: Namespace
    """
    parser = ArgumentParser(description="Split domains into train and test set based on taxonomy")

    parser.add_argument("-db", "--database", type=str, required=True,
                        help="Path to PARASECT database")
    parser.add_argument("-o", "--output", required=True, type=str,
                        help="Output directory")
    parser.add_argument("-c", '--cutoff', default=6, type=int,
                        help="Minimum substrate count for inclusion")
    parser.add_argument("-r", '--test_ratio', default=0.25, type=float,
                        help="Target test set size")
    parser.add_argument("-m", '--max_attempts', default=200, type=int,
                        help="Max attempts made for splitting data")
    parser.add_argument("-t", "--taxonomic_rank", default="family", type=str,
                        help="Taxonomic rank to split on")

    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument('--first_only', action='store_true',
                            help="Count only the first substrate per domain")
    mode_group.add_argument('--first_valid', action='store_true',
                            help="Count only the first valid substrate per domain")

    domain_group = parser.add_mutually_exclusive_group()
    domain_group.add_argument('--fungal_only', action='store_true',
                              help="Only consider fungal domains")
    domain_group.add_argument('--bacterial_only', action='store_true',
                              help="Only consider bacterial domains")

    arguments = parser.parse_args()

    return arguments


class SubstrateSelectionMode(Enum):
    FIRST_ONLY = 1
    FIRST_VALID = 2
    ALL = 3

    @staticmethod
    def get_mode_from_string(string) -> Optional["SubstrateSelectionMode"]:
        """
        Convert a string to a SplitMode enum value.

        :param string: String representation of split mode
        :type string: str
        :return: Corresponding SplitMode, or None if invalid
        :rtype: SubstrateSelectionMode
        """
        string_to_split_mode = {s.name: s for s in SubstrateSelectionMode}
        return string_to_split_mode.get(string.upper())


def _count_first_valid_substrate(session: Session, cutoff: int,
                                 included_domains: DomainScope = DomainScope.ALL) -> dict[Substrate, int]:
    """
    Count the first valid substrate for each domain in the database, iteratively filtering
    substrates that do not meet a minimum occurrence cutoff.

    This function performs the following steps:

    1. Computes the total counts of all substrates across all domains, regardless of order.
    2. Identifies a set of candidate substrates that meet the specified `cutoff` globally.
    3. Iteratively assigns to each domain the first substrate (in the domain's list) that
       is in the set of candidate substrates.
    4. After each iteration, any substrate whose assigned count falls below the `cutoff`
       is removed from the candidate set, and the process repeats until no substrates are
       pruned.

    This ensures that the returned counts only include substrates that are both globally
    and locally valid as "first substrates" for a sufficient number of domains.

    :param session: database session
    :type session: Session
    :param cutoff: minimum occurrence cutoff
    :type cutoff: int
    :param included_domains: determines which domains are considered (fungal, bacterial, all)
    :type included_domains: DomainScope
    :return: mapping of substrate to substrate counts
    :rtype: dict[Substrate, int]
    """

    # Step 1: start with all substrates meeting cutoff globally
    domains = DomainScope.get_domains(session, included_domains)
    global_counts = _count_substrates_in_domains(domains, first_only=False)
    valid_substrates = {s for s, c in global_counts.items() if c >= cutoff}

    changed = True
    counts = {}

    while changed:
        # Step 2: assign first valid substrate for each domain
        domain_to_valid = {}
        counts = {s: 0 for s in valid_substrates}
        for d in domains:
            for s in d.substrates:
                if s in valid_substrates:
                    domain_to_valid[d] = s
                    counts[s] += 1
                    break
            else:
                # No valid substrate for this domain
                domain_to_valid[d] = None

        # Step 3: prune those that fell below cutoff
        demoted = {s for s, c in counts.items() if c < cutoff}
        if demoted:
            valid_substrates -= demoted
            changed = True
        else:
            changed = False

    return counts


def _map_domains_to_first_valid(domains: Iterable[AdenylationDomain],
                                included_substrates: set[Substrate]):
    """
    Assign to each domain the first substrate from `included_substrates` that appears
    in its substrate list. Domains with no matching substrate are assigned None.

    :param domains: Domains to process.
    :type domains: Iterable[AdenylationDomain]
    :param included_substrates: Substrates considered valid for assignment.
    :type included_substrates: set[Substrate]
    :return: Mapping of domains to their first valid substrate or None.
    :rtype: dict[AdenylationDomain, Optional[Substrate]]
    """
    domain_to_first_valid = {}

    for d in domains:
        for s in d.substrates:
            if s in included_substrates:
                domain_to_first_valid[d] = s
                break
        else:
            # No valid substrate found
            domain_to_first_valid[d] = None

    return domain_to_first_valid


def get_clades(session: Session, rank: Rank,
               included_domains: DomainScope = DomainScope.ALL) -> dict[tuple[str, ...], set[AdenylationDomain]]:
    """Group domains by clades at the specified taxonomic rank, accounting for overlapping clade memberships.

    :param session: database session
    :type session: Session
    :param rank: taxonomic rank
    :type rank: Rank
    :param included_domains: determines which domains are considered (fungal, bacterial, all)
    :type included_domains: DomainScope
    :return: dictionary mapping clade groups (e.g. (Burkholderiaceae, Pseudomonadaceae) to domains
    :rtype: dict[tuple[str, ...], set[AdenylationDomain]]
    """
    # Clade grouping is necessary because it is possible that a single domain sequence maps to multiple proteins that
    # belong to different taxonomic lineages. E.g., domain ADH01485.1.A1|AIC32693.1.A1 belongs to both the
    # Burkholderiaceae (ADH01485.1.A1) and the Pseudomonadaceae (AIC32693.1.A1) families

    domains = DomainScope.get_domains(session, included_domains)
    clade_to_domains: dict[str, set[AdenylationDomain]] = {}
    clade_to_clade_group: dict[str, set[str]] = {}
    for domain in domains:
        clades = set()
        for protein in [p.protein for p in domain.proteins]:
            clades.add(Rank.get_rank_from_taxonomy(rank, protein.taxonomy))

        for clade in clades:

            # Set aliasing ensures updates across all clades
            if clade in clade_to_clade_group:
                clade_to_clade_group[clade].update(clades)
                clades = clade_to_clade_group[clade]

            else:
                clade_to_clade_group[clade] = clades

            clade_to_domains.setdefault(clade, set()).add(domain)

    # Collapse into groups: tuple of clades → set of domains
    clade_group_to_domains: dict[tuple[str, ...], set[AdenylationDomain]] = {}

    for clade, domains_in_clade in clade_to_domains.items():
        group = clade_to_clade_group[clade]
        group_key = tuple(sorted(group))
        clade_group_to_domains.setdefault(group_key, set()).update(domains_in_clade)

    return clade_group_to_domains


def _count_substrates_in_domains(domains: Iterable[AdenylationDomain],
                                 first_only: bool) -> dict[Substrate, int]:
    """Count substrates in a collection of domains.

    :param domains: collection of unique adenylation domains
    :type domains: Iterable[AdenylationDomain]
    :param first_only: if True, only count the first substrate for each domain.
    :type first_only: bool
    :return: dictionary mapping substrate to substrate counts
    :rtype: dict[Substrate, int]
    """
    counts: dict[Substrate, int] = {}

    for domain in domains:
        for substrate in domain.substrates:
            if substrate not in counts:
                counts[substrate] = 0
            counts[substrate] += 1

            if first_only:
                break

    return counts


def count_substrates(session: Session, first_only: bool, cutoff: int = 0,
                     included_domains: DomainScope = DomainScope.ALL) -> dict[Substrate, int]:
    """Count occurrences of substrates across all domains, optionally filtering by cutoff.

    :param session: database session
    :type session: Session
    :param first_only: if given, only count the first substrate for each domain
    :type first_only: bool
    :param cutoff: only include substrates that occur at least this many times
    :type cutoff: int
    :param included_domains: determines which domains are considered (fungal, bacterial, all)
    :type included_domains: DomainScope
    :return: dictionary mapping substrate to substrate counts
    :rtype: dict[Substrate, int]
    """
    domains = DomainScope.get_domains(session, included_domains)

    substrate_to_count = _count_substrates_in_domains(domains, first_only)

    for substrate, count in list(substrate_to_count.items()):
        if count < cutoff:
            del substrate_to_count[substrate]

    return substrate_to_count


def get_counts_per_clade_first_valid(clade_group_to_domains: dict[tuple[str, ...], set[AdenylationDomain]],
                                     included_substrates: set[Substrate]) -> dict[tuple[str, ...], dict[Substrate, int]]:
    """
    Count first valid substrates for each clade group.

    :param clade_group_to_domains: Mapping of clade tuples → domains
    :type clade_group_to_domains: dict[tuple[str, ...], set[AdenylationDomain]]
    :param included_substrates: Substrates considered valid
    :type included_substrates: set[Substrate]
    :return: Mapping of clade tuples → substrate counts
    :rtype: dict[tuple[str, ...], dict[Substrate, int]]
    """
    counts_per_clade = {}
    for clade, domains in clade_group_to_domains.items():
        domain_to_valid = _map_domains_to_first_valid(domains, included_substrates)
        substrate_counts = {s: 0 for s in included_substrates}
        for s in domain_to_valid.values():
            if s is not None:
                substrate_counts[s] += 1
        counts_per_clade[clade] = substrate_counts
    return counts_per_clade


def get_counts_per_clade(clade_group_to_domains: dict[tuple[str, ...], set[AdenylationDomain]],
                         included_substrates: set[Substrate],
                         first_only: bool) -> dict[tuple[str, ...], dict[Substrate, int]]:
    """Count substrates for each clade group, optionally only counting first substrates.

    :param clade_group_to_domains: dictionary mapping clade groups to domains
    :type clade_group_to_domains: dict[tuple[str, ...], set[AdenylationDomain]]
    :param included_substrates: set of included substrates
    :type included_substrates: set[Substrate]
    :param first_only: if given, only record the first substrate in the list
    :type first_only: bool
    :return: dictionary recording substrate counts for each clade
    :rtype: dict[tuple[str, ...], dict[Substrate, int]]
    """
    counts_per_clade: dict[tuple[str, ...], dict[Substrate, int]] = {}

    for clade, domains in clade_group_to_domains.items():
        substrate_to_count = _count_substrates_in_domains(domains, first_only=first_only)
        for substrate in list(substrate_to_count.keys()):
            if substrate not in included_substrates:
                del substrate_to_count[substrate]

        counts_per_clade[clade] = substrate_to_count

    return counts_per_clade


def check_split_possible(clade_substrate_counts: dict[tuple[str, ...], dict], substrates: set[Substrate]) -> bool:
    """
    Check if all substrates appear in at least two clades, required for splitting train/test sets.

    :param clade_substrate_counts: dict mapping clade groups → substrate counts
    :type clade_substrate_counts: dict[tuple[str, ...], dict]
    :param substrates: set of substrates to check
    :type substrates: set[Substrate]
    :return: True if split possible, False otherwise
    """
    for sub in substrates:
        clades_with_sub = [clade for clade, counts in clade_substrate_counts.items() if counts.get(sub, 0) > 0]
        if len(clades_with_sub) < 2:
            print(f"Substrate {sub} occurs in less than 2 clades: {len(clades_with_sub)}")
            return False
    return True


def split_on_rank(session: Session,
                  rank: Rank,
                  cutoff: int = 6,
                  target_ratio: float = 0.75,
                  split_mode: SubstrateSelectionMode = SubstrateSelectionMode.FIRST_ONLY,
                  included_domains: DomainScope = DomainScope.ALL,
                  max_attempts: int = 200) -> dict[str, dict]:
    """Split domains into train/test sets based on taxonomic rank while balancing substrate representation.

    Ensures:
      - Each substrate occurs in both train and test sets
      - Domains from the same clade are assigned to the same set
      - Substrate counts are as balanced as possible according to target_ratio

    :param session: database session
    :type session: Session
    :param rank: taxonomic rank to split on (e.g. if 'family', domains of the same family will be assigned
        to either train or test but not both
    :type rank: Rank
    :param cutoff: only consider substrates that occur at least this many times in the dataset
    :type cutoff: int
    :param target_ratio: target ratio of datapoints in training set
    :type target_ratio: float
    :param split_mode: determines which substrates are considered for splitting the data
    :type split_mode: SplitMode
    :param included_domains: determines which domains are considered (fungal, bacterial, all)
    :type split_mode: DomainScope
    :param max_attempts: maximum number of attempts to try splitting the dataset. Default: 200
    :type max_attempts: int
    :return: dictionary containing clade mappings (assignment) and substrate counts
    :rtype: dict[str, dict]
    """
    rng = random.Random(100125)
    # Decide included substrates based on mode
    if split_mode == SubstrateSelectionMode.ALL:
        substrate_to_count = count_substrates(session, first_only=False, cutoff=cutoff,
                                              included_domains=included_domains)
    elif split_mode == SubstrateSelectionMode.FIRST_ONLY:
        substrate_to_count = count_substrates(session, first_only=True, cutoff=cutoff,
                                              included_domains=included_domains)
    elif split_mode == SubstrateSelectionMode.FIRST_VALID:

        substrate_to_count = _count_first_valid_substrate(session, cutoff, included_domains=included_domains)
    else:
        raise ValueError(f"Unknown split mode: {split_mode}")

    substrates = list(substrate_to_count.keys())
    substrates.sort(key=lambda x: x.name)

    grouped_clades = get_clades(session, rank, included_domains=included_domains)

    if split_mode == SubstrateSelectionMode.FIRST_VALID:
        clade_substrate_counts = get_counts_per_clade_first_valid(grouped_clades, set(substrates))
    else:
        clade_substrate_counts = get_counts_per_clade(grouped_clades, set(substrate_to_count.keys()),
                                                      first_only=(split_mode == SubstrateSelectionMode.FIRST_ONLY))

    # coverage phase: ensure each substrate occurs in both train and test set
    def try_build_coverage():
        assignment = {}
        for substrate in rng.sample(substrates, len(substrates)):
            clades_with_substrate = [c for c, counts in clade_substrate_counts.items() if counts.get(substrate, 0) > 0]
            assigned_train = any(assignment.get(c) == "train" for c in clades_with_substrate)
            assigned_test = any(assignment.get(c) == "test" for c in clades_with_substrate)

            if assigned_train and assigned_test:
                continue

            if not assigned_train and not assigned_test:
                if len(clades_with_substrate) < 2:
                    return None
                c1, c2 = rng.sample(clades_with_substrate, 2)
                assignment[c1] = "train"
                assignment[c2] = "test"

            elif assigned_train and not assigned_test:
                unassigned = [c for c in clades_with_substrate if c not in assignment]
                if not unassigned:
                    return None
                assignment[unassigned[0]] = "test"

            elif assigned_test and not assigned_train:
                unassigned = [c for c in clades_with_substrate if c not in assignment]
                if not unassigned:
                    return None
                assignment[unassigned[0]] = "train"

        return assignment

    coverage_assignment = None

    for _ in range(max_attempts):
        coverage_assignment = try_build_coverage()
        if coverage_assignment:
            break
    if not coverage_assignment:
        raise RuntimeError("Could not satisfy coverage constraints.")

    # initialize substrate counts
    substrate_counts = {s: {"train": 0, "test": 0} for s in substrates}
    for clade, side in coverage_assignment.items():
        for sub, cnt in clade_substrate_counts[clade].items():
            substrate_counts[sub][side] += cnt

    # assign remaining clades greedily
    unassigned = [c for c in clade_substrate_counts if c not in coverage_assignment]
    unassigned.sort(key=lambda c: -sum(clade_substrate_counts[c].values()))

    def total_error(counts):
        err = 0.0
        for sub, ct in counts.items():
            total = ct["train"] + ct["test"]
            if total == 0:
                continue
            desired = total * target_ratio
            err += (ct["train"] - desired) ** 2
        return err

    for clade in unassigned:
        best_side, best_err = None, None
        for side in ("train", "test"):
            temp = {s: ct.copy() for s, ct in substrate_counts.items()}
            for sub, cnt in clade_substrate_counts[clade].items():
                temp[sub][side] += cnt
            err = total_error(temp)
            if best_err is None or err < best_err:
                best_side, best_err = side, err
        coverage_assignment[clade] = best_side
        for sub, cnt in clade_substrate_counts[clade].items():
            substrate_counts[sub][best_side] += cnt

    # sanity check
    for sub, ct in substrate_counts.items():
        if ct["train"] == 0 or ct["test"] == 0:
            raise RuntimeError(f"Substrate {sub} missing from one split")

    domain_assignment: dict[AdenylationDomain, str] = {}
    for clade, side in coverage_assignment.items():
        for domain in grouped_clades[clade]:
            if domain in domain_assignment:
                raise RuntimeError(f"Domain {domain} assigned twice")

            if split_mode == SubstrateSelectionMode.ALL:
                if any(sub in substrates for sub in domain.substrates):
                    domain_assignment[domain] = side
            elif split_mode == SubstrateSelectionMode.FIRST_ONLY:
                if domain.substrates[0] in substrates:
                    domain_assignment[domain] = side
            elif split_mode == SubstrateSelectionMode.FIRST_VALID:
                # assign based on first valid substrate
                for s in domain.substrates:
                    if s in substrates:
                        domain_assignment[domain] = side
                        break

    return {
        "assignment": coverage_assignment,
        "substrate_counts": substrate_counts,
        "domain_assignment": domain_assignment
    }


def main():
    args = parse_arguments()

    if args.first_only:
        split_mode = SubstrateSelectionMode.FIRST_ONLY
    elif args.first_valid:
        split_mode = SubstrateSelectionMode.FIRST_VALID
    else:
        split_mode = SubstrateSelectionMode.ALL

    if args.fungal_only:
        included_domains = DomainScope.FUNGAL_ONLY
    elif args.bacterial_only:
        included_domains = DomainScope.BACTERIAL_ONLY
    else:
        included_domains = DomainScope.ALL

    engine = create_engine(f"sqlite:///{args.database}")
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    train_clades_out = os.path.join(args.output, "train_clades.txt")
    test_clades_out = os.path.join(args.output, "test_clades.txt")

    train_out = os.path.join(args.output, "train.txt")
    test_out = os.path.join(args.output, "test.txt")

    substrate_splits = os.path.join(args.output, "substrate_splits.txt")

    with Session(engine) as session:
        rank_type = Rank.get_rank_type_from_string(args.taxonomic_rank)
        result = split_on_rank(session, rank_type,
                               cutoff=args.cutoff,
                               target_ratio=1.0 - args.test_ratio,
                               split_mode=split_mode,
                               included_domains=included_domains,
                               max_attempts=args.max_attempts)

        clade_assignment = result["assignment"]
        substrate_counts = result["substrate_counts"]
        domain_assignment = result["domain_assignment"]

        with open(train_out, 'w') as train, open(test_out, 'w') as test:
            for domain, assignment in domain_assignment.items():
                if assignment == "train":
                    train.write(f"{domain.get_name()}\n")
                elif assignment == "test":
                    test.write(f"{domain.get_name()}\n")
                else:
                    raise ValueError(f"Unknown domain assignment: {assignment}")

        with open(train_clades_out, 'w') as train_clades, open(test_clades_out, 'w') as test_clades:
            for clades, assignment in clade_assignment.items():
                for clade in clades:
                    if assignment == "train":
                        train_clades.write(f"{clade}\n")
                    elif assignment == "test":
                        test_clades.write(f"{clade}\n")
                    else:
                        raise ValueError(f"Unknown clade assignment: {assignment}")

        with open(substrate_splits, 'w') as substrates_out:
            substrates_out.write("substrate_name\ttrain\ttest\n")
            for substrate, counts in substrate_counts.items():
                substrates_out.write(f"{substrate.name}\t{counts['train']}\t{counts['test']}\n")

        # Sum across all substrates
        total_train = sum(ct["train"] for ct in substrate_counts.values())
        total_test = sum(ct["test"] for ct in substrate_counts.values())

        print(f"Total domains in train: {total_train}")
        print(f"Total domains in test: {total_test}")

        if args.first_only:
            num_train_domains = sum(1 for side in domain_assignment.values() if side == "train")
            num_test_domains = sum(1 for side in domain_assignment.values() if side == "test")

            total_counts = total_train + total_test
            num_domains = len(domain_assignment)

            # Check totals
            if total_counts != num_domains:
                raise RuntimeError(
                    f"Sanity check failed: {total_counts} substrate counts but "
                    f"{num_domains} domains assigned (expected equality with first_only=True)."
                )

            # Check per-split counts
            if total_train != num_train_domains:
                raise RuntimeError(
                    f"Sanity check failed: train substrate count ({total_train}) "
                    f"!= number of train domains ({num_train_domains})."
                )

            if total_test != num_test_domains:
                raise RuntimeError(
                    f"Sanity check failed: test substrate count ({total_test}) "
                    f"!= number of test domains ({num_test_domains})."
                )
        elif not args.first_valid:
            # 1. No domain assigned twice
            assigned_domains = set()
            for domain, side in domain_assignment.items():
                if domain in assigned_domains:
                    raise RuntimeError(f"Domain {domain} assigned twice")
                assigned_domains.add(domain)

            # 2. Every substrate counted is on the side of some domain
            for sub, counts in substrate_counts.items():
                counted_train = sum(1 for d, side in domain_assignment.items()
                                    if side == "train" and sub in d.substrates)
                counted_test = sum(1 for d, side in domain_assignment.items()
                                   if side == "test" and sub in d.substrates)

                if counted_train != counts["train"]:
                    raise RuntimeError(
                        f"Mismatch for substrate {sub}: train counted {counts['train']}, "
                        f"sum over domains {counted_train}"
                    )
                if counted_test != counts["test"]:
                    raise RuntimeError(
                        f"Mismatch for substrate {sub}: test counted {counts['test']}, "
                        f"sum over domains {counted_test}"
                    )

            # Get all domains that have at least one substrate in the counted list
            legal_domains = set()
            for domain in DomainScope.get_domains(session, included_domains):
                if any(sub in substrate_counts for sub in domain.substrates):
                    legal_domains.add(domain)

            # Compare with assigned domains
            assigned_domains = set(domain_assignment.keys())

            missing = legal_domains - assigned_domains
            extra = assigned_domains - legal_domains

            if missing or extra:
                raise RuntimeError(
                    f"Domain assignment mismatch: "
                    f"{len(missing)} legal domains missing, "
                    f"{len(extra)} extra domains assigned"
                )

        else:
            assigned_domains = set()

            # 1. Ensure no domain assigned twice
            for domain, side in domain_assignment.items():
                if domain in assigned_domains:
                    raise RuntimeError(f"Domain {domain} assigned twice")
                assigned_domains.add(domain)

            # 2. Verify counted substrates match assigned domains
            for sub, counts in substrate_counts.items():
                counted_train = sum(1 for d, side in domain_assignment.items()
                                    if
                                    side == "train" and sub == next((s for s in d.substrates if s in substrate_counts),
                                                                    None))
                counted_test = sum(1 for d, side in domain_assignment.items()
                                   if side == "test" and sub == next((s for s in d.substrates if s in substrate_counts),
                                                                     None))

                if counted_train != counts["train"]:
                    raise RuntimeError(
                        f"Mismatch for substrate {sub}: train counted {counts['train']}, sum over domains {counted_train}"
                    )
                if counted_test != counts["test"]:
                    raise RuntimeError(
                        f"Mismatch for substrate {sub}: test counted {counts['test']}, sum over domains {counted_test}"
                    )

            # 3. Ensure all domains with at least one valid substrate are assigned
            legal_domains = set()
            for domain in DomainScope.get_domains(session, included_domains):
                if any(sub in substrate_counts for sub in domain.substrates):
                    legal_domains.add(domain)

            missing = legal_domains - assigned_domains
            extra = assigned_domains - legal_domains

            if missing or extra:
                raise RuntimeError(
                    f"Domain assignment mismatch: "
                    f"{len(missing)} legal domains missing, "
                    f"{len(extra)} extra domains assigned"
                )


if __name__ == "__main__":
    main()
