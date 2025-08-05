from typing import Dict, Type, Union

from mat3ra.made.material import Material
from .configuration import NanoparticleConfiguration
from ... import MaterialWithBuildMetadata
from ...defect.point.vacancy.builder import VacancyDefectBuilder
from ...merge import MergeBuilder
from ...slab.slab.builder import SlabBuilder
from ...slab.slab.configuration import SlabConfiguration
from ...void_region.builder import VoidRegionBuilder
from ...void_region.configuration import VoidRegionConfiguration
from ....analyze.other import get_chemical_formula_empirical


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
