from statistics import mean

from paras.scripts.feature_extraction.sequence_feature_extraction.seq_to_features import PROPERTIES


def get_gap_properties():
    index_to_properties: dict[int, list[float]] = {}
    for aa, property_vector in PROPERTIES.items():
        if aa != '-':
            for i, value in enumerate(property_vector):
                if i not in index_to_properties:
                    index_to_properties[i] = []
                index_to_properties[i].append(value)

    index_to_average: dict[int, float] = {}
    for index, properties in index_to_properties.items():
        assert len(properties) == 20
        index_to_average[index] = mean(properties)

    return index_to_average


if __name__ == "__main__":
    index_to_average = get_gap_properties()
    property_string = '-'
    for i in range(len(index_to_average.keys())):
        property_string += f'\t{index_to_average[i]:.3f}'

    print(property_string)
