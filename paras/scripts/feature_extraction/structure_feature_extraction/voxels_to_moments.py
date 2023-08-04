from dataclasses import asdict
from collections import OrderedDict
from paras.scripts.feature_extraction.structure_feature_extraction.moment_invariants import get_moments
from paras.scripts.parsers.voxel import voxels_from_file

from sys import argv

MOMENT_TYPES = ['O_3', 'O_4', 'O_5', 'F', 'phi_2', 'phi_3', 'phi_4', 'phi_5', 'phi_6', 'phi_7', 'phi_8',
                'phi_9', 'phi_10', 'phi_11', 'phi_12', 'phi_13', 'CI']


class MomentDataset:

    def __init__(self, voxel_file):
        self.domain_to_moments = OrderedDict()
        self.attributes = []
        attributes_set = False
        counter = 0
        for domain_name, voxels in voxels_from_file(voxel_file):
            attribute_to_moments = self.voxels_to_moments(voxels)
            if not attributes_set:
                self.attributes = list(attribute_to_moments.keys())
                attributes_set = True

            moment_vector = []
            for attribute, moments in attribute_to_moments.items():
                moment_vector += moments

            self.domain_to_moments[domain_name] = moment_vector
            counter += 1
            if counter % 10 == 0:
                print(f"Calculated moments for {counter} pockets.")

    @staticmethod
    def voxels_to_moments(voxels):
        attribute_to_coords = OrderedDict()
        for voxel in voxels:
            attribute_to_counts = asdict(voxel.vector)
            for attribute, count in attribute_to_counts.items():
                if attribute not in attribute_to_coords:
                    attribute_to_coords[attribute] = []
                if count:
                    attribute_to_coords[attribute].append(voxel.cube.midpoint.to_list())

        attributes = sorted(attribute_to_coords.keys())
        attribute_to_moments = OrderedDict()
        for attribute in attributes:
            coordinates = attribute_to_coords[attribute]
            if len(coordinates) < 2:
                moment_invariants = [0.0] * 17
            else:
                moment_invariants = get_moments(coordinates)
            attribute_to_moments[attribute] = list(moment_invariants)

        return attribute_to_moments

    def write_moments(self, out_file):
        with open(out_file, 'w') as out:
            out.write("domain_name\t")
            column_labels = []
            for attribute in self.attributes:
                for moment_type in MOMENT_TYPES:
                    column_label = f"{attribute}_{moment_type}"
                    column_labels.append(column_label)

            out.write('\t'.join(column_labels))
            out.write('\n')

            for domain, moments in self.domain_to_moments.items():
                out.write(f"{domain}")
                for moment in moments:
                    out.write(f"\t{moment:.10f}")
                out.write('\n')


if __name__ == "__main__":
    moment_data = MomentDataset(argv[1])
    moment_data.write_moments(argv[2])
