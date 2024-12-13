from enum import Enum


class NanoparticleShapes(str, Enum):
    """
    Enum for supported nanoparticle shapes.
    """

    CUBOCTAHEDRON = "cuboctahedron"
    ICOSAHEDRON = "icosahedron"
    OCTAHEDRON = "octahedron"
    DECAHEDRON = "decahedron"
    SIMPLE_CUBIC = "simple_cubic"
    BODY_CENTERED_CUBIC = "body_centered_cubic"
    HEXAGONAL_CLOSED_PACKED = "hexagonal_closed_packed"
    WULFF = "wulff"
