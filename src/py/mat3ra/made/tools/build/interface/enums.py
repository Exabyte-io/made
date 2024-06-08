from enum import Enum


class StrainModes(str, Enum):
    strain = "strain"
    von_mises_strain = "von_mises_strain"
    mean_abs_strain = "mean_abs_strain"
