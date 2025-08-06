from typing import Dict, Type, Union

from mat3ra.made.material import Material
from .configuration import NanoparticleConfiguration
from .....defective_structures.zero_dimensional.point_defect.vacancy.builder import VacancyDefectBuilder
from ......analyze import get_chemical_formula_empirical
from ......build_components import MaterialWithBuildMetadata
from ......build_components.entities.auxiliary.zero_dimensional.void_region.builder import VoidRegionBuilder
from ......build_components.entities.auxiliary.zero_dimensional.void_region.configuration import VoidRegionConfiguration
from mat3ra.made.tools.build.pristine_structures.two_dimensional.slab.builder import SlabBuilder
from mat3ra.made.tools.build.pristine_structures.two_dimensional.slab.configuration import SlabConfiguration
from ......build_components.operations.core.combinations.merge import MergeBuilder


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
