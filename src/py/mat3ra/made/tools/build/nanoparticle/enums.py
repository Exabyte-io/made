from enum import Enum
from typing import Callable
from ...third_party import (
    ASEAtoms,
    ASESimpleCubic,
    ASEBodyCenteredCubic,
    ASEFaceCenteredCubic,
    ASEIcosahedron,
    ASEOctahedron,
    ASEDecahedron,
    ASEHexagonalClosedPacked,
    ASEWulffConstruction,
)


class ASENanoparticleShapesEnum(str, Enum):
    """
    Enum for supported nanoparticle shapes.
    """

    ICOSAHEDRON = "icosahedron"
    OCTAHEDRON = "octahedron"
    DECAHEDRON = "decahedron"
    SIMPLE_CUBIC = "simple_cubic"
    FACE_CENTERED_CUBIC = "face_centered_cubic"
    BODY_CENTERED_CUBIC = "body_centered_cubic"
    HEXAGONAL_CLOSED_PACKED = "hexagonal_closed_packed"
    WULFF_CONSTRUCTION = "wulff_construction"

    @staticmethod
    def get_ase_constructor(shape: str) -> Callable[..., ASEAtoms]:
        return {
            "icosahedron": ASEIcosahedron,
            "octahedron": ASEOctahedron,
            "decahedron": ASEDecahedron,
            "simple_cubic": ASESimpleCubic,
            "face_centered_cubic": ASEFaceCenteredCubic,
            "body_centered_cubic": ASEBodyCenteredCubic,
            "hexagonal_closed_packed": ASEHexagonalClosedPacked,
            "ase_wulff_construction": ASEWulffConstruction,
        }[shape]
