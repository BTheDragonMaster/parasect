class DataPoint:
    def __init__(self):
        self.id = ''
        self.sequence = ''
        self.specificities = set()
        self.major_specificities = set()
        self.minor_specificities = set()
        self.specificity_to_activity = {}
        self.evidence = set()
        self.source = ''
        
    def set_source(self, source: str) -> None:
        self.source = source

    def set_sequence(self, sequence: str) -> None:
        self.sequence = sequence

    def set_id(self, seq_id: str) -> None:
        self.id = seq_id

    def add_specificity(self, specificity: str) -> None:
        self.specificities.add(specificity)
        self.major_specificities.add(specificity)

    def add_minor_specificity(self, specificity: str) -> None:
        self.minor_specificities.add(specificity)

    def add_evidence(self, evidence: str) -> None:
        self.evidence.add(evidence)

    def set_activity(self, specificity: str, activity: float) -> None:
        assert 0.0 <= activity <= 1.0
        if specificity not in self.specificity_to_activity:
            self.specificity_to_activity[specificity] = activity