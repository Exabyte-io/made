from typing import List, Callable, Dict, Type, Union

from mat3ra.made.material import Material
from .configuration import ASEBasedNanoparticleConfiguration, NanoparticleConfiguration
from .enums import NanoparticleShapesEnum
from .. import MaterialWithBuildMetadata
from ..defect.point.builders import VacancyDefectBuilder
from ..merge import MergeBuilder
from ..slab.builders import SlabBuilder
from ..slab.configurations import SlabConfiguration
from ..void_region.builders import VoidRegionBuilder
from ..void_region.configuration import VoidRegionConfiguration
from ...analyze.other import get_chemical_formula_empirical
from ...build import BaseBuilder
from ...build.mixins import ConvertGeneratedItemsASEAtomsMixin
from ...third_party import ASEAtoms


class NanoparticleBuilder(VacancyDefectBuilder, MergeBuilder):
    """
    Builder class for creating a nanoparticle by merging a slab with a void region.
    """

    @property
    def merge_component_types_conversion_map(self) -> Dict[Type, Type]:
        return {
            SlabConfiguration: SlabBuilder,
            VoidRegionConfiguration: VoidRegionBuilder,
        }

    def _update_material_name(
        self, material: Union[Material, MaterialWithBuildMetadata], configuration: NanoparticleConfiguration
    ) -> Material:
        formula = get_chemical_formula_empirical(material)

        material.name = f"{formula} {configuration.void_region_configuration.condition_name} Nanoparticle"
        return material


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
        return parameters

    @classmethod
    def _get_ase_nanoparticle_constructor(cls, config: ASEBasedNanoparticleConfiguration) -> Callable[..., ASEAtoms]:
        constructor = NanoparticleShapesEnum.get_ase_constructor(config.shape.value)
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
            material.name = f"{get_chemical_formula_empirical(material)} {configuration.shape.value.capitalize()}"
        return materials
