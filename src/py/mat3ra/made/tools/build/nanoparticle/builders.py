from typing import Callable, Dict, Type, Union

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
from ...build import BaseSingleBuilder
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


class ASEBasedNanoparticleBuilder(ConvertGeneratedItemsASEAtomsMixin, BaseSingleBuilder):
    """
    Generalized builder for creating nanoparticles based on ASE cluster tools.
    Passes configuration parameters directly to the ASE constructors.
    """

    _ConfigurationType: type(ASEBasedNanoparticleConfiguration) = ASEBasedNanoparticleConfiguration  # type: ignore
    _GeneratedItemType: type(ASEAtoms) = ASEAtoms  # type: ignore

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
    def _set_ase_cell(atoms: ASEAtoms, vacuum: float) -> ASEAtoms:
        """
        Set the cell of an ASE atoms object to a cubic box with vacuum padding around the nanoparticle.
        """
        max_dimension_along_x = abs(atoms.positions[:, 0]).max()
        box_size = 2 * max_dimension_along_x + vacuum
        atoms.set_cell([box_size, box_size, box_size], scale_atoms=False)
        atoms.center()
        return atoms

    def _generate(self, configuration: ASEBasedNanoparticleConfiguration) -> MaterialWithBuildMetadata:
        parameters = self._get_ase_nanoparticle_parameters(configuration)
        constructor = self._get_ase_nanoparticle_constructor(configuration)
        nanoparticle_without_cell = constructor(**parameters)
        nanoparticle = self._set_ase_cell(nanoparticle_without_cell, configuration.vacuum_padding)

        return MaterialWithBuildMetadata.create(self._convert_generated_item(nanoparticle))

    def _update_material_name(
        self, material: MaterialWithBuildMetadata, configuration: ASEBasedNanoparticleConfiguration
    ) -> MaterialWithBuildMetadata:
        material.name = f"{get_chemical_formula_empirical(material)} {configuration.shape.value.capitalize()}"
        return material
