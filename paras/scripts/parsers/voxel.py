from dataclasses import dataclass, astuple
from typing import List, Optional, Dict
from paras.scripts.math.shapes import Cube, Vector3D
from paras.scripts.parsers.atom import Atom
from paras.scripts.parsers.tabular import Tabular
from enum import Enum, unique
from paras.scripts.data_processing.structure_processing.voxel_adjacency_matrix import make_adjacency_matrix

from sys import argv


ATOM_TO_RESIDUE_TO_TYPE = {'CG': {"Y": ['aromatic'],
                                  "F": ['aromatic'],
                                  "W": ['aromatic'],
                                  "H": ['aromatic']},
                           'CD1': {"Y": ['aromatic'],
                                   "F": ['aromatic'],
                                   "W": ['aromatic']},
                           'CD2': {"Y": ['aromatic'],
                                   "F": ['aromatic'],
                                   "W": ['aromatic'],
                                   "H": ['aromatic']},
                           'CE1': {"Y": ['aromatic'],
                                   "F": ['aromatic'],
                                   "W": ['aromatic'],
                                   "H": ['aromatic']},
                           'CE2': {"Y": ['aromatic'],
                                   "F": ['aromatic'],
                                   "W": ['aromatic']},
                           'CE3': {"W": ['aromatic']},
                           'CZ': {"Y": ['aromatic'],
                                  "F": ['aromatic'],
                                  "W": ['aromatic']},
                           'CZ2': {"W": ['aromatic']},
                           'CZ3': {"W": ['aromatic']},
                           'CH2': {"W": ['aromatic']},
                           'ND1': {"W": ['aromatic'],
                                   "H": ['aromatic', '+']},
                           'ND2': {"N": ['pos_polar']},
                           'NE1': {"W": ['aromatic', 'pos_polar']},
                           'NE2': {"Q": ['pos_polar'],
                                   "H": ['aromatic', '+']},
                           'OD1': {"D": ['-'],
                                   'N': ['neg_polar']},
                           'OD2': {"D": ['-']},
                           'OE1': {"E": ['-'],
                                   "Q": ['neg_polar']},
                           'OE2': {"E": ['-']},
                           'NE': {"R": ['+']},
                           'NH1': {"R": ['+']},
                           'NH2': {"R": ['+']},
                           'NZ': {"K": ['+']},
                           'OG': {"S": ['neg_polar']},
                           'OG1': {"T": ['neg_polar']},
                           'SD': {"M": ['sulphur']},
                           'SG': {'C': ['sulphur']},
                           'OH': {'Y': ['neg_polar']},
                           'P': {'AMP': ['phosphorus']},
                           'O1P': {'AMP': ['-']},
                           'O2P': {'AMP': ['-']},
                           'O3P': {'AMP': ['-']},
                           "O5'": {'AMP': ['neg_polar']},
                           "O4'": {'AMP': ['neg_polar']},
                           "O3'": {'AMP': ['neg_polar']},
                           "O2'": {'AMP': ['neg_polar']},
                           "N1": {'AMP': ['pos_polar', 'aromatic']},
                           "N3": {'AMP': ['pos_polar', 'aromatic']},
                           "N7": {'AMP': ['pos_polar', 'aromatic']},
                           "N9": {'AMP': ['pos_polar', 'aromatic']},
                           "N6": {'AMP': ['pos_polar']},
                           "C2": {'AMP': ['aromatic']},
                           "C4": {'AMP': ['aromatic']},
                           "C5": {'AMP': ['aromatic']},
                           "C6": {'AMP': ['aromatic']},
                           "C8": {'AMP': ['aromatic']},
                           "MG": {'MG': ['metal']}
                           }


VOXEL_TYPE_TO_ABBREVIATION = {"hydrophobic": "APOL",
                              "oxygen": "NEGA",
                              "sulphur": "CYST",
                              "phosphorus": "PHOS",
                              "nitrogen": "POSI",
                              "polar": " POL",
                              "metal": " MET",
                              "aromatic": "AROM"}


@unique
class VoxelDataType(Enum):
    is_empty = 1
    is_pocket = 2
    nr_aromatic = 3
    nr_hydrophobic = 4
    nr_pos_charge = 5
    nr_neg_charge = 6
    nr_pos_polar = 7
    nr_neg_polar = 8
    nr_sulphur = 9
    nr_phosphorus = 10
    nr_metal = 11

    @staticmethod
    def from_string(label: str) -> "VoxelDataType":
        for value in VoxelDataType:
            if str(value.name) == label:
                return value
        raise ValueError(f"Unknown voxel data type: {label}")


@dataclass
class VoxelVector:
    is_empty: int = 0
    is_pocket: int = 0
    nr_aromatic: int = 0
    nr_hydrophobic: int = 0
    nr_pos_charge: int = 0
    nr_neg_charge: int = 0
    nr_pos_polar: int = 0
    nr_neg_polar: int = 0
    nr_sulfur: int = 0
    nr_phosphorus: int = 0
    nr_metal: int = 0

    def update_vector(self, category_name, category_value):
        setattr(self, category_name, category_value)

    def as_list(self, categories: Optional[List] = None, data_type: str = "int"):
        if categories is None:
            attribute_list = [getattr(self, "is_empty"),
                              getattr(self, "is_pocket"),
                              getattr(self, "nr_aromatic"),
                              getattr(self, "nr_hydrophobic"),
                              getattr(self, "nr_pos_charge"),
                              getattr(self, "nr_neg_charge"),
                              getattr(self, "nr_pos_polar"),
                              getattr(self, "nr_neg_polar"),
                              getattr(self, "nr_sulfur"),
                              getattr(self, "nr_phosphorus"),
                              getattr(self, "nr_metal")]
        else:
            attribute_list = []
            for category in categories:
                attribute_list.append(getattr(self, category))

        if data_type == 'float':
            return list(map(float, attribute_list))
        elif data_type == 'int':
            return attribute_list
        else:
            raise ValueError(f"Expected string 'int' or 'float'. Got {data_type}.")


@dataclass
class VoxelImportanceVector:
    is_empty: float = 0.0
    is_pocket: float = 0.0
    nr_aromatic: float = 0.0
    nr_hydrophobic: float = 0.0
    nr_pos_charge: float = 0.0
    nr_neg_charge: float = 0.0
    nr_pos_polar: float = 0.0
    nr_neg_polar: float = 0.0
    nr_sulfur: float = 0.0
    nr_phosphorus: float = 0.0
    nr_metal: float = 0.0

    def update_vector(self, category_name, category_value):
        setattr(self, category_name, category_value)


@dataclass
class Voxel:
    cube: Cube
    atoms: List[Atom]
    type: Optional[str] = None
    polar: bool = False
    pocket: bool = False
    vector: Optional[VoxelVector] = None
    importances: Optional[VoxelImportanceVector] = None

    def __post_init__(self):
        if self.vector is None:
            self.vector = VoxelVector()

    def __hash__(self):
        return hash(self.cube.midpoint)

    def __eq__(self, other):
        if self.cube.midpoint == other.cube.midpoint:
            return True
        else:
            return False

    def __repr__(self):
        return self.cube.midpoint.__repr__()

    def update_vector(self, category_name, category_value):
        self.vector.update_vector(category_name, category_value)

    def update_importances(self, category_name, category_value):
        self.importances.update_vector(category_name, category_value)

    def to_feature_vector(self):

        atom_type_to_count = {"aromatic": 0,
                              "-": 0,
                              '+': 0,
                              'pos_polar': 0,
                              'neg_polar': 0,
                              'sulphur': 0,
                              'hydrophobic': 0,
                              'metal': 0,
                              'phosphorus': 0}
        for atom in self.atoms:
            if atom.name in ATOM_TO_RESIDUE_TO_TYPE:
                residue_to_type = ATOM_TO_RESIDUE_TO_TYPE[atom.name]
                if atom.residue.type in residue_to_type:
                    for atom_type in residue_to_type[atom.residue.type]:
                        atom_type_to_count[atom_type] += 1
                elif atom.type == 'C':
                    atom_type_to_count['hydrophobic'] += 1
                else:
                    print(atom.name)
            elif atom.type == 'C':
                atom_type_to_count['hydrophobic'] += 1
            elif atom.name == 'N':
                atom_type_to_count['pos_polar'] += 1
            elif atom.name in {'O', 'OXT'}:
                atom_type_to_count['neg_polar'] += 1
            else:
                print(atom.name)

        if self.type == 'empty':
            empty = 1
        else:
            empty = 0

        if self.pocket:
            pocket = 1
        else:
            pocket = 0

        self.vector = VoxelVector(empty,
                                  pocket,
                                  atom_type_to_count["aromatic"],
                                  atom_type_to_count["+"],
                                  atom_type_to_count["-"],
                                  atom_type_to_count["pos_polar"],
                                  atom_type_to_count["neg_polar"],
                                  atom_type_to_count["sulphur"],
                                  atom_type_to_count["hydrophobic"],

                                  atom_type_to_count["phosphorus"],
                                  atom_type_to_count["metal"])

        return list(astuple(self.vector))


def voxels_from_file(voxel_file, voxel_size=1.2, pocket_only=True):
    with open(voxel_file, 'r') as voxel_data:
        categories = voxel_data.readline().strip().split('\t')
        for line in voxel_data:
            voxel_to_vector = {}
            voxel_info = line.strip().split('\t')
            domain_name = voxel_info[0]
            for i, voxel_value in enumerate(voxel_info):
                if i == 0:
                    continue
                category = categories[i]
                coord_string, category_name = category.split('|')
                coords = list(map(float, coord_string.split('_')))
                cube = Cube(Vector3D(coords[0], coords[1], coords[2]), voxel_size)
                voxel = Voxel(cube, [])
                if voxel not in voxel_to_vector:
                    voxel_to_vector[voxel] = VoxelVector()
                voxel_to_vector[voxel].update_vector(category_name, int(voxel_value))

            for voxel, vector in voxel_to_vector.items():
                voxel.vector = vector

            voxels = list(voxel_to_vector.keys())

            if pocket_only:
                pocket_voxels = []
                for voxel in voxels:
                    if voxel.vector.is_pocket:
                        pocket_voxels.append(voxel)
                adjacency_matrix = make_adjacency_matrix(pocket_voxels, voxels, voxel_size)
                pocket_voxels = set()
                for voxel, neighbours in adjacency_matrix.items():
                    pocket_voxels.add(voxel)
                    for neighbour in neighbours:
                        pocket_voxels.add(neighbour)

                for voxel in pocket_voxels:
                    voxel.update_vector('is_pocket', 1)
                yield domain_name, pocket_voxels
            else:
                yield domain_name, voxels


def change_headers(voxel_file, voxels_out):
    with open(voxels_out, 'w') as out:
        with open(voxel_file, 'r') as voxels:
            header = voxels.readline()
            header = header.replace('i#sulphur', 'nr_sulfur')
            header = header.replace('phosphorus', 'nr_phosphorus')
            header = header.replace('metal', 'nr_metal')
            header = header.replace('hydrophobic', 'nr_hydrophobic')
            header = header.replace('#+', 'nr_pos_charge')
            header = header.replace('#-', 'nr_neg_charge')
            header = header.replace('#pos_polar', 'nr_pos_polar')
            header = header.replace('#neg_polar', 'nr_neg_polar')
            header = header.replace('#aromatic', 'nr_aromatic')
            out.write(header)
            for line in voxels:
                out.write(line)


if __name__ == "__main__":
    for domain_name, voxels in voxels_from_file(argv[1], pocket_only=True):
        print(domain_name, [voxel.vector for voxel in voxels])
        print(len(voxels))
        exit()
