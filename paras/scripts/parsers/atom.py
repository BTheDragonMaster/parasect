from dataclasses import dataclass
from paras.scripts.math.shapes import Sphere


@dataclass
class Residue:
    type: str
    nr: int

    def __hash__(self):
        return hash((self.type, self.nr))

    def __eq__(self, other):
        if self.type == other.type and self.nr == other.nr:
            return True
        return False


@dataclass
class Atom:
    type: str
    name: str
    residue: Residue
    chain_id: str
    sphere: Sphere

    def __repr__(self):
        return f"{self.name}_{self.residue.type}{self.residue.nr}"
