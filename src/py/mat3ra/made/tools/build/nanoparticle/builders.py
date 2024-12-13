from typing import List
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
from mat3ra.made.material import Material
from mat3ra.made.tools.build import BaseBuilder
from mat3ra.made.tools.convert import from_ase

from .configuration import NanoparticleConfiguration
from .enums import NanoparticleShapes


class NanoparticleBuilder(BaseBuilder):
    """
    Generalized builder for creating nanoparticles based on ASE cluster tools.
    Passes configuration parameters directly to the ASE constructors.
    """

    _ConfigurationType: type(NanoparticleConfiguration) = NanoparticleConfiguration  # type: ignore
    _GeneratedItemType: Material = Material

    def create_nanoparticle(self, config: NanoparticleConfiguration) -> Material:
        shape = config.shape
        element = config.element

        lattice_constant = config.lattice_constant

        # Ensure parameters dictionary exists
        parameters = config.parameters or {}

        # Add common parameters
        parameters["symbol"] = element
        if "latticeconstant" not in parameters:
            parameters["latticeconstant"] = lattice_constant

        # Shape-specific factory logic
        if shape == NanoparticleShapes.CUBOCTAHEDRON:
            nanoparticle = FaceCenteredCubic(**parameters)
        elif shape == NanoparticleShapes.ICOSAHEDRON:
            nanoparticle = Icosahedron(**parameters)
        elif shape == NanoparticleShapes.OCTAHEDRON:
            nanoparticle = Octahedron(**parameters)
        elif shape == NanoparticleShapes.DECAHEDRON:
            nanoparticle = Decahedron(**parameters)
        elif shape == NanoparticleShapes.SIMPLE_CUBIC:
            nanoparticle = SimpleCubic(**parameters)
        elif shape == NanoparticleShapes.BODY_CENTERED_CUBIC:
            nanoparticle = BodyCenteredCubic(**parameters)
        elif shape == NanoparticleShapes.HEXAGONAL_CLOSED_PACKED:
            nanoparticle = HexagonalClosedPacked(**parameters)
        elif shape == NanoparticleShapes.WULFF:
            nanoparticle = wulff_construction(**parameters)
        else:
            raise ValueError(f"Unsupported shape: {shape}")

        return Material(from_ase(nanoparticle))

    def _generate(self, configuration: NanoparticleConfiguration) -> List[_GeneratedItemType]:
        nanoparticle = self.create_nanoparticle(configuration)
        return [nanoparticle]
