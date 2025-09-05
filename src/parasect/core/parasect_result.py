from typing import Union

from parasect.core.domain import AdenylationDomain
from parasect.core.constants import SEPARATOR_1, SEPARATOR_2, SEPARATOR_3


class Result:
    """Result class."""

    def __init__(
        self,
        domain: AdenylationDomain,
        predictions: list[float],
        prediction_labels: list[str],
        prediction_smiles: list[str],
    ) -> None:
        """Initialise the Result class.

        :param domain: Adenylation domain.
        :type domain: AdenylationDomain
        :param predictions: Predictions.
        :type predictions: List[float]
        :param prediction_labels: Prediction labels.
        :type prediction_labels: List[str]
        :param prediction_smiles: Prediction SMILES.
        :type prediction_smiles: List[str]
        """
        self.domain = domain
        self._predictions = predictions
        self._prediction_labels = prediction_labels
        self._prediction_smiles = prediction_smiles

    @property
    def predictions(self):
        return self._predictions

    @property
    def prediction_labels(self):
        return self._prediction_labels

    def sort(self):
        pred_label_smiles = list(zip(self._predictions, self._prediction_labels, self._prediction_smiles))
        pred_label_smiles.sort(key=lambda x: x[0], reverse=True)
        self._predictions = [data[0] for data in pred_label_smiles]
        self._prediction_labels = [data[1] for data in pred_label_smiles]
        self._prediction_smiles = [data[2] for data in pred_label_smiles]

    def get_domain_header(self, separator_1: str = SEPARATOR_1, separator_2: str = SEPARATOR_2,
                          separator_3: str = SEPARATOR_3) -> str:
        """Get domain header for domain

        :param separator_1: Separator between header fields
        :type separator_1: str
        :param separator_2: Separator between 'domain' string and domain number
        :type separator_2: str
        :param separator_3: Separator between domain start and end coordinates
        :type separator_3: str
        :return: domain header
        :rtype: str
        """
        domain_header = f"{self.domain.protein_name}{separator_1}domain{separator_2}{self.domain.domain_nr}{separator_1}{self.domain.start}{separator_3}{self.domain.end}"
        return domain_header

    def to_json(self) -> dict[str, Union[str, int, list[dict[str, Union[str, float]]]]]:
        """Return the Result as a JSON serialisable dictionary.

        :return: JSON serialisable dictionary.
        :rtype: Dict[str, Union[str, List[(float, str)]]
        """
        return dict(
            domain_name=self.domain.protein_name,
            domain_nr=self.domain.domain_nr,
            domain_start=self.domain.start,
            domain_end=self.domain.end,
            domain_sequence=self.domain.sequence,
            domain_signature=self.domain.signature,
            domain_extended_signature=self.domain.extended_signature,
            predictions=[
                dict(substrate_name=sub_name, substrate_smiles=sub_smiles, probability=prob)
                for sub_name, sub_smiles, prob in zip(
                    self._prediction_labels,
                    self._prediction_smiles,
                    self._predictions,
                )
            ],
        )
