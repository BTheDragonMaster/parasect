SEPARATOR_1 = '|'
SEPARATOR_2 = '_'
SEPARATOR_3 = '-'


def parse_domain_id(fasta_id):
    """
    Return id, domain type and domain location from common id

    Input:
    fasta_id: str

    Output:
    id: str, sequence id
    hit_id: str, domain id
    hit_start: int, start position of domain in protein
    hit_end: int, end position of domain in protei
    """
    seq_id, hit_id, hit_location = fasta_id.split(SEPARATOR_1)
    hit_start, hit_end = hit_location.split('-')
    hit_start = int(hit_start)
    hit_end = int(hit_end)
    return seq_id, hit_id, hit_start, hit_end


def make_domain_id(seq_id, hit_id, start, end):
    """
    Return header (str) for .fasta output file

    Input:
    ID: str, sequence identifier
    hit_id: str, HMM the sequence matched to
    start: int, first aa in the query sequence that the HMM matches to
    end: int, last aa in the query sequence that the HMM matches to

    Output:
    str, sequence properties separated by '|'
    """
    return f'{seq_id}{SEPARATOR_1}{hit_id}{SEPARATOR_1}{start}-{end}'
