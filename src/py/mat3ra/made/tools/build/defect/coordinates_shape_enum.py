from enum import Enum


class CoordinatesShapeEnum(str, Enum):
    SPHERE = "sphere"
    CYLINDER = "cylinder"
    BOX = "rectangle"
    TRIANGULAR_PRISM = "triangular_prism"
