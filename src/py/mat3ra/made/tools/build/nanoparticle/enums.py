from enum import Enum
from typing import Callable
from ...third_party import ASEAtoms

# TODO: import from 3rd party
from ase.cluster import (
    SimpleCubic,
    BodyCenteredCubic,
    FaceCenteredCubic,
    Icosahedron,
    Octahedron,
    Decahedron,
    HexagonalClosedPacked,
)
from ase.cluster.wulff import wulff_construction


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
            "icosahedron": Icosahedron,
            "octahedron": Octahedron,
            "decahedron": Decahedron,
            "simple_cubic": SimpleCubic,
            "face_centered_cubic": FaceCenteredCubic,
            "body_centered_cubic": BodyCenteredCubic,
            "hexagonal_closed_packed": HexagonalClosedPacked,
            "wulff_construction": wulff_construction,
        }[shape]
