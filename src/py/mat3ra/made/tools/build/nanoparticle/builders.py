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
from mat3ra.made.tools.build.mixins import ConvertGeneratedItemsASEAtomsMixin
from mat3ra.made.tools.build.slab import SlabConfiguration
from mat3ra.made.tools.modify import filter_by_condition_on_coordinates

from .configuration import ASENanoparticleConfiguration, NanoparticleConfiguration
from .enums import NanoparticleShapes
from ..slab import create_slab
from ...third_party import ASEAtoms


class NanoparticleBuilder(ConvertGeneratedItemsASEAtomsMixin, BaseBuilder):
    """
    Builder for creating nanoparticles by cutting from bulk materials supercells.
    """

    _ConfigurationType: type(NanoparticleConfiguration) = NanoparticleConfiguration  # type: ignore
    _GeneratedItemType: type(Material) = Material  # type: ignore

    def create_nanoparticle(self, config: _ConfigurationType) -> _GeneratedItemType:
        material = config.material
        # shape = config.shape
        orientation_z = config.orientation_z
        radius = config.radius

        # Get the conventional structure
        conventional_material = self._ConfigurationType.convert_to_conventional(material)
        thickness_in_layers = (2 * radius + config.vacuum_padding) / conventional_material.lattice_constant
        slab_config = SlabConfiguration(
            bulk=conventional_material.structure,
            miller_indices=orientation_z,
            thickness=thickness_in_layers,
            use_conventional_cell=True,
            use_orthogonal_z=True,
            make_primitive=False,
        )
        slab = create_slab(slab_config)
        # TODO: switch condition based on shape
        nanoparticle = filter_by_condition_on_coordinates(slab, lambda vector: vector.norm() <= radius)
        return nanoparticle


class ASEBasedNanoparticleBuilder(ConvertGeneratedItemsASEAtomsMixin, NanoparticleBuilder):
    """
    Generalized builder for creating nanoparticles based on ASE cluster tools.
    Passes configuration parameters directly to the ASE constructors.
    """

    _ConfigurationType: type(ASENanoparticleConfiguration) = ASENanoparticleConfiguration  # type: ignore
    _GeneratedItemType: type(ASEAtoms) = ASEAtoms  # type: ignore

    def create_nanoparticle(self, config: ASENanoparticleConfiguration) -> _GeneratedItemType:
        shape = config.shape
        element = config.element

        lattice_constant = config.lattice_constant

        # Ensure parameters dictionary exists
        parameters = config.parameters or {}

        # Add common parameters
        parameters["symbol"] = element
        if "latticeconstant" not in parameters:
            parameters["latticeconstant"] = lattice_constant
        # TODO: adjust parameters for octahedron based on type (cuboctahedron, regular octahedron, etc.)
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

    def _generate(self, configuration: ASENanoparticleConfiguration) -> List[_GeneratedItemType]:
        nanoparticle = self.create_nanoparticle(configuration)
        return [nanoparticle]
