from paras.scripts.feature_extraction.sequence_feature_extraction.read_positions import \
    HMM2_POSITIONS_EXTENDED_SIGNATURE, HMM2_POSITIONS_SIGNATURE, HMM2_POSITION_K, POSITIONS_SIGNATURE, \
    POSITIONS_EXTENDED_SIGNATURE, \
    get_reference_positions_hmm, get_reference_positions

from paras.scripts.feature_extraction.sequence_feature_extraction.profile_alignment.align import align_adomain, \
    ALIGNMENT_FILE

VALID_CHARACTERS = {"A", "C", "D", "E",
                    "F", "G", "H", "I",
                    "K", "L", "M", "N",
                    "P", "Q", "R", "S",
                    "T", "V", "W", "Y",
                    "-"}


class AdenylationDomain:
    def __init__(self, protein_name, domain_start, domain_end):
        self.protein_name = protein_name
        self.start = domain_start
        self.end = domain_end

        self.sequence = None
        self.signature = None
        self.extended_signature = None
        self.domain_nr = 0
        self.domain_id = None

    def set_domain_id(self, separator_1, separator_2, separator_3):
        assert self.protein_name
        assert self.domain_nr

        self.domain_id = f"{self.protein_name}{separator_1}domain{separator_2}{self.domain_nr}{separator_1}{self.start}{separator_3}{self.end}"

    def set_domain_signatures_hmm(self, hit_n_terminal, hit_c_terminal=None):
        """Extract (extended) signatures from adenylation domains using HMM"""

        signature_positions = HMM2_POSITIONS_SIGNATURE
        extended_signature_positions = HMM2_POSITIONS_EXTENDED_SIGNATURE
        position_k = HMM2_POSITION_K

        profile = hit_n_terminal.aln[1].seq
        query = hit_n_terminal.aln[0].seq
        offset = hit_n_terminal.hit_start

        signature = get_reference_positions_hmm(query, profile, [p - offset for p in signature_positions])
        if signature and all([char in VALID_CHARACTERS for char in signature]):
            self.signature = signature

        lysine = None

        if hit_c_terminal:
            profile_c = hit_c_terminal.aln[1].seq
            query_c = hit_c_terminal.aln[0].seq
            offset_c = hit_c_terminal.hit_start
            lysine = get_reference_positions_hmm(query_c, profile_c, [p - offset_c for p in position_k])

        if self.signature:
            if lysine and lysine in VALID_CHARACTERS and lysine != '-':
                self.signature += lysine
            else:
                self.signature += "K"

        extended_signature = get_reference_positions_hmm(query, profile, [p - offset for p in extended_signature_positions])
        if extended_signature and all([char in VALID_CHARACTERS for char in extended_signature]):
            self.extended_signature = extended_signature

    def set_domain_signatures_profile(self):
        """Extract (extended) signatures from adenylation domains using profile alignment"""

        if not self.sequence:
            raise Exception("Sequence needs to be defined first")

        seq_id = "DOMAIN_TO_QUERY"

        aligned_domain, aligned_reference = align_adomain(seq_id, self.sequence, ALIGNMENT_FILE)

        aligned_positions_signature = get_reference_positions(POSITIONS_SIGNATURE, aligned_reference)
        aligned_positions_extended_signature = get_reference_positions(POSITIONS_EXTENDED_SIGNATURE, aligned_reference)

        signature = []
        for position in aligned_positions_signature:
            signature.append(aligned_domain[position])
        self.signature = ''.join(signature)

        extended_signature = []
        for position in aligned_positions_extended_signature:
            extended_signature.append(aligned_domain[position])
        self.extended_signature = ''.join(extended_signature)

    def set_domain_number(self, domain_nr):
        self.domain_nr = domain_nr

    def set_sequence(self, sequence):
        self.sequence = sequence

    def write_sequence(self, fasta_file, sequence_type='full'):
        assert self.domain_id
        if sequence_type == 'full':
            fasta_file.write(f'>{self.domain_id}\n{self.sequence}\n')
        elif sequence_type == 'extended_signature':
            fasta_file.write(f'>{self.domain_id}\n{self.extended_signature}\n')
        elif sequence_type == 'signature':
            fasta_file.write(f'>{self.domain_id}\n{self.signature}\n')
        else:
            raise ValueError(f"Only accepted sequence types are 'full', 'signature', and 'extended_signature'. Got {sequence_type}")
