import os.path
from argparse import ArgumentParser, Namespace

from parasect.core.featurisation import get_domains
from parasect.core.parsing import parse_fasta_file, write_fasta_file


def parse_arguments() -> Namespace:
    """
    Parse arguments for extracting A domain signatures

    :return: arguments
    :rtype: Namespace
    """
    parser = ArgumentParser(description="Extract signatures for A-domain dataset")
    parser.add_argument('-i', required=True, type=str, help="Path to input .fasta")
    parser.add_argument('-t', required=True, type=str, help="Path to temp dir")
    parser.add_argument('-o', required=True, type=str, help="Path to output directory")

    args = parser.parse_args()

    return args


def extract_signatures(fasta_in: str, tmp_path: str, out_path: str) -> None:
    """Extract signatures for A-domain dataset

    :param fasta_in: path to fasta input file containing adenylation domains
    :type fasta_in: str
    :param tmp_path: path to tmp parasect path for storing intermediate files
    :type tmp_path: str
    :param out_path: path to output directory
    :type out_path: str

    """
    original_domains = set(parse_fasta_file(fasta_in).keys())
    domains = get_domains(fasta_in, tmp_path, "hmm", "fasta")

    if not os.path.exists(out_path):
        os.mkdir(out_path)

    signatures_out = os.path.join(out_path, 'signatures.fasta')
    extended_out = os.path.join(out_path, 'extended_signatures.fasta')

    for domain in domains:
        if domain.domain_nr != 1:
            print(domain.protein_name)

    id_to_sig: dict[str, str] = {}
    id_to_ext: dict[str, str] = {}

    for domain in domains:
        original_domain = domain.protein_name
        if original_domain in id_to_sig:
            print(f"Warning: multiple domains found in {original_domain}. Signatures: {id_to_sig[original_domain]}, {domain.signature}. Selecting signature with fewest gaps. If this is an A-OX domain, revise signature manually.")
            if id_to_sig[original_domain].count('-') > domain.signature.count('-'):
                id_to_sig[original_domain] = domain.signature
                id_to_ext[original_domain] = domain.extended_signature
        else:
            id_to_sig[original_domain] = domain.signature
            id_to_ext[original_domain] = domain.extended_signature

        if original_domain in original_domains:
            original_domains.remove(original_domain)

    for domain in original_domains:
        print(f"Could not extract signatures for {domain}. Attempt extraction with full protein.")

    write_fasta_file(id_to_sig, signatures_out)
    write_fasta_file(id_to_ext, extended_out)


def main() -> None:
    args = parse_arguments()
    extract_signatures(args.i, args.t, args.o)


if __name__ == "__main__":
    main()
