from typing import List, Callable, Dict, Type

from mat3ra.made.material import Material
from ...analyze.other import get_chemical_formula
from ...build import BaseBuilder
from ...build.mixins import ConvertGeneratedItemsASEAtomsMixin
from ...build.slab import SlabConfiguration
from ...modify import filter_by_condition_on_coordinates
from ...utils.coordinate import SphereCoordinateCondition
from ...analyze.other import get_closest_site_id_from_coordinate
from ...third_party import ASEAtoms
from ..slab import create_slab
from .configuration import ASEBasedNanoparticleConfiguration, SphereSlabBasedNanoparticleConfiguration
from .enums import ASENanoparticleShapesEnum


class SlabBasedNanoparticleBuilder(BaseBuilder):
    """
    Builder for creating nanoparticles by cutting from bulk materials supercells.
    """

    _ConfigurationType: type(ASEBasedNanoparticleConfiguration) = ASEBasedNanoparticleConfiguration  # type: ignore
    _GeneratedItemType: type(Material) = Material  # type: ignore

    def create_nanoparticle(self, config: _ConfigurationType) -> _GeneratedItemType:
        slab = self._create_slab(config)
        center_coordinate = self._find_slab_center_coordinate(slab)
        condition = self._build_coordinate_condition(config, center_coordinate)
        nanoparticle = filter_by_condition_on_coordinates(slab, condition, use_cartesian_coordinates=True)
        return nanoparticle

    def _build_coordinate_condition(self, config: _ConfigurationType, center_coordinate: List[float]) -> Callable:
        coordinate_condition = config.condition_builder(center_coordinate)
        return coordinate_condition.condition

    def _create_slab(self, config: _ConfigurationType) -> Material:
        slab_config = SlabConfiguration(
            bulk=config.material,
            miller_indices=config.orientation_z,
            thickness=config.supercell_size,
            use_conventional_cell=True,
            use_orthogonal_z=True,
            make_primitive=False,
            vacuum=0,
            xy_supercell_matrix=[[config.supercell_size, 0], [0, config.supercell_size]],
        )
        slab = create_slab(slab_config)
        return slab

    def _find_slab_center_coordinate(self, slab: Material) -> List[float]:
        slab.to_cartesian()
        center_coordinate = slab.basis.cell.convert_point_to_cartesian([0.5, 0.5, 0.5])
        center_id_at_site = get_closest_site_id_from_coordinate(slab, center_coordinate, use_cartesian_coordinates=True)
        center_coordinate_at_site = slab.basis.coordinates.get_element_value_by_index(center_id_at_site)
        return center_coordinate_at_site

    def _finalize(self, materials: List[Material], configuration: _ConfigurationType) -> List[Material]:
        for material in materials:
            material.name = f"{get_chemical_formula(material)} Nanoparticle"
        return materials


class SphereSlabBasedNanoparticleBuilder(SlabBasedNanoparticleBuilder):
    """
    Builder for creating spherical nanoparticles by cutting from bulk materials supercells.
    """

    _ConfigurationType: Type[ASEBasedNanoparticleConfiguration] = SphereSlabBasedNanoparticleConfiguration

    def _build_coordinate_condition(
        self, config: SphereSlabBasedNanoparticleConfiguration, center_coordinate: List[float]
    ) -> Callable:
        return SphereCoordinateCondition(center_position=center_coordinate, radius=config.radius).condition


class ASEBasedNanoparticleBuilder(ConvertGeneratedItemsASEAtomsMixin, BaseBuilder):
    """
    Generalized builder for creating nanoparticles based on ASE cluster tools.
    Passes configuration parameters directly to the ASE constructors.
    """

    _ConfigurationType: type(ASEBasedNanoparticleConfiguration) = ASEBasedNanoparticleConfiguration  # type: ignore
    _GeneratedItemType: type(ASEAtoms) = ASEAtoms  # type: ignore

    def create_nanoparticle(self, config: ASEBasedNanoparticleConfiguration) -> _GeneratedItemType:
        parameters = self._get_ase_nanoparticle_parameters(config)
        constructor = self._get_ase_nanoparticle_constructor(config)
        nanoparticle_without_cell = constructor(**parameters)
        nanoparticle = self._set_ase_cell(nanoparticle_without_cell, config.vacuum_padding)

        return nanoparticle

    @staticmethod
    def _get_ase_nanoparticle_parameters(config: ASEBasedNanoparticleConfiguration) -> Dict:
        parameters = config.parameters or {}
        parameters["symbol"] = config.element
        parameters.setdefault("latticeconstant", config.lattice_constant)
        return parameters

    @classmethod
    def _get_ase_nanoparticle_constructor(cls, config: ASEBasedNanoparticleConfiguration) -> Callable[..., ASEAtoms]:
        constructor = ASENanoparticleShapesEnum.get_ase_constructor(config.shape.value)
        return constructor

    @staticmethod
    def _set_ase_cell(atoms: ASEAtoms, vacuum: float) -> _GeneratedItemType:
        """
        Set the cell of an ASE atoms object to a cubic box with vacuum padding around the nanoparticle.
        """
        max_dimension_along_x = abs(atoms.positions[:, 0]).max()
        box_size = 2 * max_dimension_along_x + vacuum
        atoms.set_cell([box_size, box_size, box_size], scale_atoms=False)
        atoms.center()
        return atoms

    def _generate(self, configuration: ASEBasedNanoparticleConfiguration) -> List[_GeneratedItemType]:
        nanoparticle = self.create_nanoparticle(configuration)
        return [nanoparticle]

    def _finalize(self, materials: List[Material], configuration: _ConfigurationType) -> List[Material]:
        for material in materials:
            material.name = f"{get_chemical_formula(material)} {configuration.shape.value.capitalize()}"
        return materials
