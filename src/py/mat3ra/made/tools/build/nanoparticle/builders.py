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
from mat3ra.made.tools.build import BaseBuilder
from mat3ra.made.tools.build.mixins import ConvertGeneratedItemsASEAtomsMixin

from .configuration import NanoparticleConfiguration
from .enums import NanoparticleShapes
from ...third_party import ASEAtoms


class NanoparticleBuilder(ConvertGeneratedItemsASEAtomsMixin, BaseBuilder):
    """
    Generalized builder for creating nanoparticles based on ASE cluster tools.
    Passes configuration parameters directly to the ASE constructors.
    """

    _ConfigurationType: type(NanoparticleConfiguration) = NanoparticleConfiguration  # type: ignore
    _GeneratedItemType: type(ASEAtoms) = ASEAtoms  # type: ignore

    def create_nanoparticle(self, config: NanoparticleConfiguration) -> _GeneratedItemType:
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
        box_size = 2 * max(abs(nanoparticle.positions).max(axis=0)) + config.vacuum_padding
        nanoparticle.set_cell([box_size, box_size, box_size], scale_atoms=False)
        nanoparticle.center()
        return nanoparticle

    def _generate(self, configuration: NanoparticleConfiguration) -> List[_GeneratedItemType]:
        nanoparticle = self.create_nanoparticle(configuration)
        return [nanoparticle]
