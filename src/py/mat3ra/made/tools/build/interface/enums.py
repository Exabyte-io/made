from enum import Enum
from typing import TypedDict, List


class StrainModes(str, Enum):
    strain = "strain"
    von_mises_strain = "von_mises_strain"
    mean_abs_strain = "mean_abs_strain"


class SupercellTypes(str, Enum):
    hexagonal = "hexagonal"
    orthogonal = "orthogonal"


class SupercellMatrix(TypedDict):
    angle: float
    xy_supercell: List[List[int]]


# Tabulation from https://github.com/qtm-iisc/Twister/blob/474156a2a59f2b9d59350b32de56864a9496f848/examples/Homobilayer_hex/hex.table
# Maps twist angle to supercell matrix values for homo-material hexagonal supercells bilayer
angle_to_supercell_matrix_values_for_hex: List[SupercellMatrix] = [
    {"angle": 60.0, "xy_supercell": [[0, 1], [-1, 1]]},
    {"angle": 21.7867892983, "xy_supercell": [[1, 2], [-2, 3]]},
    {"angle": 13.1735511073, "xy_supercell": [[2, 3], [-3, 5]]},
    {"angle": 9.4300079079, "xy_supercell": [[3, 4], [-4, 7]]},
    {"angle": 7.34099301663, "xy_supercell": [[4, 5], [-5, 9]]},
    {"angle": 6.00898319777, "xy_supercell": [[5, 6], [-6, 11]]},
    {"angle": 5.08584780812, "xy_supercell": [[6, 7], [-7, 13]]},
    {"angle": 4.40845500794, "xy_supercell": [[7, 8], [-8, 15]]},
    {"angle": 3.89023816901, "xy_supercell": [[8, 9], [-9, 17]]},
    {"angle": 3.48100608947, "xy_supercell": [[9, 10], [-10, 19]]},
    {"angle": 3.14965742639, "xy_supercell": [[10, 11], [-11, 21]]},
    {"angle": 2.87589463363, "xy_supercell": [[11, 12], [-12, 23]]},
    {"angle": 2.64590838119, "xy_supercell": [[12, 13], [-13, 25]]},
    {"angle": 2.44997727662, "xy_supercell": [[13, 14], [-14, 27]]},
    {"angle": 2.28105960973, "xy_supercell": [[14, 15], [-15, 29]]},
    {"angle": 2.1339296666, "xy_supercell": [[15, 16], [-16, 31]]},
    {"angle": 2.00462783069, "xy_supercell": [[16, 17], [-17, 33]]},
    {"angle": 1.89009907347, "xy_supercell": [[17, 18], [-18, 35]]},
    {"angle": 1.78794861038, "xy_supercell": [[18, 19], [-19, 37]]},
    {"angle": 1.69627269352, "xy_supercell": [[19, 20], [-20, 39]]},
    {"angle": 1.61353890116, "xy_supercell": [[20, 21], [-21, 41]]},
    {"angle": 1.53849981837, "xy_supercell": [[21, 22], [-22, 43]]},
    {"angle": 1.47012972578, "xy_supercell": [[22, 23], [-23, 45]]},
    {"angle": 1.40757744635, "xy_supercell": [[23, 24], [-24, 47]]},
    {"angle": 1.35013073557, "xy_supercell": [[24, 25], [-25, 49]]},
    {"angle": 1.29718904759, "xy_supercell": [[25, 26], [-26, 51]]},
    {"angle": 1.24824246553, "xy_supercell": [[26, 27], [-27, 53]]},
    {"angle": 1.20285522748, "xy_supercell": [[27, 28], [-28, 55]]},
    {"angle": 1.16065271985, "xy_supercell": [[28, 29], [-29, 57]]},
    {"angle": 1.12131111538, "xy_supercell": [[29, 30], [-30, 59]]},
    {"angle": 1.08454904916, "xy_supercell": [[30, 31], [-31, 61]]},
    {"angle": 1.05012087979, "xy_supercell": [[31, 32], [-32, 63]]},
    {"angle": 1.01781119445, "xy_supercell": [[32, 33], [-33, 65]]},
    {"angle": 0.987430297814, "xy_supercell": [[33, 34], [-34, 67]]},
    {"angle": 0.958810485525, "xy_supercell": [[34, 35], [-35, 69]]},
    {"angle": 0.931802947264, "xy_supercell": [[35, 36], [-36, 71]]},
    {"angle": 0.906275178895, "xy_supercell": [[36, 37], [-37, 73]]},
    {"angle": 0.882108808579, "xy_supercell": [[37, 38], [-38, 75]]},
    {"angle": 0.859197761631, "xy_supercell": [[38, 39], [-39, 77]]},
    {"angle": 0.8374467041, "xy_supercell": [[39, 40], [-40, 79]]},
    {"angle": 0.816769716893, "xy_supercell": [[40, 41], [-41, 81]]},
    {"angle": 0.0, "xy_supercell": [[1, 0], [0, 1]]},
]
