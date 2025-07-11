from typing import Any

from mat3ra.made.material import Material
from mat3ra.made.tools.build.defect.slab.configuration import SlabDefectConfiguration, AdatomDefectConfiguration, \
    IslandDefectConfiguration, VoidSite
from mat3ra.made.tools.build.merge.builders import MergeBuilder
from mat3ra.made.tools.modify import filter_by_condition_on_coordinates
from ..slab.configuration import SlabStackConfiguration
from ...slab.builders import SlabBuilder
from ...slab.configurations import SlabConfiguration
from ...stack.builders import StackNComponentsBuilder


class SlabStackBuilder(StackNComponentsBuilder):
    def _configuration_to_material(self, configuration_or_material: Any) -> Material:
        if isinstance(configuration_or_material, SlabConfiguration):
            return SlabBuilder().get_material(configuration_or_material)
        else:
            return super()._configuration_to_material(configuration_or_material)

    def _generate(self, configuration: SlabStackConfiguration) -> Material:
        return super()._generate(configuration)



class IslandDefectBuilder(MergeBuilder):
    """
    Builder for creating island defects by merging a slab with a void site.
    The void site defines which atoms to remove from the slab.
    """

    _ConfigurationType = IslandDefectConfiguration

    def _configuration_to_material(self, configuration_or_material):
        if isinstance(configuration_or_material, VoidSite):
            return self._create_void_material(configuration_or_material)
        return super()._configuration_to_material(configuration_or_material)

    def _create_void_material(self, void_site: VoidSite) -> Material:
        """
        Create a material that represents the void by removing atoms that don't satisfy the condition.
        """
        if not void_site.crystal:
            raise ValueError("VoidSite must have a crystal defined")

        # Create a copy of the crystal
        void_material = void_site.crystal.clone()

        # Remove atoms that don't satisfy the condition
        filtered_material = filter_by_condition_on_coordinates(
            material=void_material,
            condition=void_site.condition_function,
            use_cartesian_coordinates=False,  # Use crystal coordinates by default
        )

        return filtered_material
