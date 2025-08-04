from enum import Enum


class BondDirectionsTemplatesEnum(list, Enum):
    OCTAHEDRAL = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]]
    TETRAHEDRAL = [[1, 1, 1], [-1, -1, 1], [1, -1, -1], [-1, 1, -1]]
    PLANAR = [[1, 0, 0], [-1, 2, 0], [-1, -2, 0]]
    LINEAR = [[1, 0, 0], [-1, 0, 0]]
