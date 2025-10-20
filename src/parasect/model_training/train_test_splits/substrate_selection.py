from typing import Optional, Iterable, Union
from enum import Enum

from sqlalchemy.orm import Session

from parasect.database.build_database import Substrate, AdenylationDomain

from .domain_scope import DomainScope


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


def count_first_valid_substrate(session: Session, cutoff: int,
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
    global_counts = count_substrates_in_domains(domains, first_only=False)
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


def count_substrates_in_domains(domains: Iterable[AdenylationDomain],
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

    substrate_to_count = count_substrates_in_domains(domains, first_only)

    for substrate, count in list(substrate_to_count.items()):
        if count < cutoff:
            del substrate_to_count[substrate]

    return substrate_to_count


def map_domains_to_substrates(domains: Iterable[AdenylationDomain],
                              included_substrates: set[Substrate],
                              selection_mode: SubstrateSelectionMode) -> \
        dict[AdenylationDomain, Union[Optional[list[Substrate]], Optional[Substrate]]]:
    """
    Assign to each domain one or more substrates. Domains with no matching substrate are assigned None.

    :param domains: Domains to process.
    :type domains: Iterable[AdenylationDomain]
    :param included_substrates: Substrates considered valid for assignment.
    :type included_substrates: set[Substrate]
    :param selection_mode: substrate selection mode (consider the first valid substrate, only the first substrate,
    or all substrates
    :return: Mapping of domains to their substrate(s) or None.
    :rtype: dict[AdenylationDomain, Union[Optional[list[Substrate]], Optional[Substrate]]]
    """

    domain_to_substrate = {}

    for domain in domains:
        if selection_mode == SubstrateSelectionMode.FIRST_VALID:
            for substrate in domain.substrates:
                if substrate in included_substrates:
                    domain_to_substrate[domain] = substrate
                    break
            else:
                # No valid substrate found
                domain_to_substrate[domain] = None
        elif selection_mode == SubstrateSelectionMode.FIRST_ONLY:
            substrate = domain.substrates[0]
            if substrate in included_substrates:
                domain_to_substrate[domain] = substrate
            else:
                domain_to_substrate[domain] = None
        elif selection_mode == SubstrateSelectionMode.ALL:
            domain_to_substrate[domain] = []
            for substrate in domain.substrates:
                if substrate in included_substrates:
                    domain_to_substrate[domain].append(substrate)

            if not domain_to_substrate[domain]:
                domain_to_substrate[domain] = None

    return domain_to_substrate
