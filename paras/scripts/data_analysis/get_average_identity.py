from statistics import mean, stdev
from sys import argv

from paras.scripts.parsers.parsers import parse_closest_identities, parse_domain_list


def get_average_identity(identity_file, domain_list):
    domain_to_closest = parse_closest_identities(identity_file)
    included_domains = parse_domain_list(domain_list)
    identities = []

    for domain in included_domains:
        identity = domain_to_closest[domain]
        identities.append(identity)

    print(f"Mean: {mean(identities)}")
    print(f"Standard deviation: {stdev(identities)}")


if __name__ == "__main__":
    get_average_identity(argv[1], argv[2])

