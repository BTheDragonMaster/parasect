import numpy as np
from sklearn.preprocessing import MultiLabelBinarizer

def binarise_data(labels: list[list[str]]) -> tuple[np.ndarray, np.ndarray]:
    mlb = MultiLabelBinarizer()
    new_labels = mlb.fit_transform(labels)
    return new_labels, mlb.classes_