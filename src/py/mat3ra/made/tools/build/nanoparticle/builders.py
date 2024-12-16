from typing import List, Callable, Dict
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
from mat3ra.made.tools.analyze.other import get_chemical_formula
from mat3ra.made.tools.build import BaseBuilder
from mat3ra.made.tools.build.mixins import ConvertGeneratedItemsASEAtomsMixin
from mat3ra.made.tools.build.slab import SlabConfiguration
from mat3ra.made.tools.modify import filter_by_condition_on_coordinates

from .configuration import ASEBasedNanoparticleConfiguration, NanoparticleConfiguration
from .enums import NanoparticleShapes
from ..slab import create_slab
from ...third_party import ASEAtoms

SHAPE_TO_CONSTRUCTOR: Dict[str, Callable[..., ASEAtoms]] = {
    NanoparticleShapes.ICOSAHEDRON: Icosahedron,
    NanoparticleShapes.OCTAHEDRON: Octahedron,
    NanoparticleShapes.DECAHEDRON: Decahedron,
    NanoparticleShapes.SIMPLE_CUBIC: SimpleCubic,
    NanoparticleShapes.FACE_CENTERED_CUBIC: FaceCenteredCubic,
    NanoparticleShapes.BODY_CENTERED_CUBIC: BodyCenteredCubic,
    NanoparticleShapes.HEXAGONAL_CLOSED_PACKED: HexagonalClosedPacked,
    NanoparticleShapes.WULFF: wulff_construction,
}


class NanoparticleBuilder(BaseBuilder):
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

    def _finalize(self, materials: List[Material], configuration: _ConfigurationType) -> List[Material]:
        for material in materials:
            material.name = f"{get_chemical_formula(material)} {configuration.shape.value.capitalize()}"
        return materials


class ASEBasedNanoparticleBuilder(ConvertGeneratedItemsASEAtomsMixin, NanoparticleBuilder):
    """
    Generalized builder for creating nanoparticles based on ASE cluster tools.
    Passes configuration parameters directly to the ASE constructors.
    """

    _ConfigurationType: type(ASEBasedNanoparticleConfiguration) = ASEBasedNanoparticleConfiguration  # type: ignore
    _GeneratedItemType: type(ASEAtoms) = ASEAtoms  # type: ignore

    def create_nanoparticle(self, config: ASEBasedNanoparticleConfiguration) -> _GeneratedItemType:
        constructor = SHAPE_TO_CONSTRUCTOR.get(config.shape.value)
        if not constructor:
            raise ValueError(f"Unsupported shape: {config.shape}")
        parameters = config.parameters or {}
        parameters["symbol"] = config.element
        parameters.setdefault("latticeconstant", config.lattice_constant)
        nanoparticle = constructor(**parameters)
        box_size = 2 * max(abs(nanoparticle.positions).max(axis=0)) + config.vacuum_padding
        nanoparticle.set_cell([box_size, box_size, box_size], scale_atoms=False)
        nanoparticle.center()

        return nanoparticle

    def _generate(self, configuration: ASEBasedNanoparticleConfiguration) -> List[_GeneratedItemType]:
        nanoparticle = self.create_nanoparticle(configuration)
        return [nanoparticle]

    def _finalize(self, materials: List[Material], configuration: _ConfigurationType) -> List[Material]:
        for material in materials:
            material.name = f"{get_chemical_formula(material)} {configuration.shape.value.capitalize()}"
        return materials
