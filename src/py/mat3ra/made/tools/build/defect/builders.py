from typing import List, Callable, Optional, Union

from mat3ra.made.material import Material
from mat3ra.made.utils import get_atomic_coordinates_extremum
from .configuration import (
    PointDefectConfigurationLegacy,
    AdatomSlabPointDefectConfiguration,
    IslandSlabDefectConfiguration,
    TerraceSlabDefectConfiguration,
    PointDefectPairConfiguration,
)
from .factories import DefectBuilderFactory
from .point.builders import PointDefectBuilder
from ...build import BaseBuilder
from ...modify import (
    filter_by_box,
    filter_by_condition_on_coordinates,
)


class DefectBuilder(BaseBuilder):
    def create_isolated_defect(self, defect_configuration: PointDefectConfigurationLegacy) -> Material:
        raise NotImplementedError


class DefectPairBuilder(DefectBuilder):
    def create_defect_pair(
        self,
        primary_defect_configuration: Union[PointDefectConfigurationLegacy, AdatomSlabPointDefectConfiguration],
        secondary_defect_configuration: Union[PointDefectConfigurationLegacy, AdatomSlabPointDefectConfiguration],
    ) -> Material:
        """
        Create a pair of point defects in the material.

        Args:
            primary_defect_configuration: The configuration of the first defect.
            secondary_defect_configuration: The configuration of the second defect.

        Returns:
            Material: The material with both defects added.
        """
        primary_material = self.create_isolated_defect(primary_defect_configuration)
        # Remove metadata to allow for independent defect creation
        if hasattr(primary_defect_configuration.crystal.metadata, "build"):
            primary_material.metadata["build"] = primary_defect_configuration.crystal.metadata["build"]
        primary_material.name = primary_defect_configuration.crystal.name
        secondary_defect_configuration.crystal = primary_material
        secondary_material = self.create_isolated_defect(secondary_defect_configuration)

        return secondary_material


class PointDefectPairBuilder(PointDefectBuilder, DefectPairBuilder):
    _ConfigurationType: type(PointDefectPairConfiguration) = PointDefectPairConfiguration  # type: ignore
    _GeneratedItemType: Material = Material

    def create_isolated_defect(self, defect_configuration: PointDefectConfigurationLegacy) -> Material:
        key = defect_configuration.defect_type.value
        if hasattr(defect_configuration, "placement_method") and defect_configuration.placement_method is not None:
            key += f":{defect_configuration.placement_method.name}".lower()
        builder_class = DefectBuilderFactory.get_class_by_name(key)
        defect_builder = builder_class()
        return defect_builder.get_material(defect_configuration)

    def _generate(self, configuration: _ConfigurationType) -> List[_GeneratedItemType]:
        return [
            self.create_defect_pair(
                primary_defect_configuration=configuration.primary_defect_configuration,
                secondary_defect_configuration=configuration.secondary_defect_configuration,
            )
        ]

    def _update_material_name(self, material: Material, configuration: _ConfigurationType) -> Material:
        updated_material = super()._update_material_name(material, configuration)
        name_1 = configuration.primary_defect_configuration.defect_type.name.capitalize()
        name_2 = configuration.secondary_defect_configuration.defect_type.name.capitalize()
        new_name = f"{updated_material.name}, {name_1} and {name_2} Defect Pair"
        updated_material.name = new_name
        return updated_material


class IslandSlabDefectBuilder(DefectBuilder):
    _ConfigurationType: type(IslandSlabDefectConfiguration) = IslandSlabDefectConfiguration  # type: ignore
    _GeneratedItemType: Material = Material

    @staticmethod
    def _default_condition(coordinate: List[float]):
        return True

    def create_island(
        self,
        material: Material,
        condition: Optional[Callable[[List[float]], bool]] = None,
        thickness: int = 1,
        use_cartesian_coordinates: bool = False,
    ) -> Material:
        """
        Create an island at the specified position on the surface of the material.

        Args:
            material: The material to add the island to.
            condition: The condition on coordinates to describe the island.
            thickness: The thickness of the island in layers.
            use_cartesian_coordinates: Whether to use Cartesian coordinates for the condition.
        Returns:
            The material with the island added.
        """
        new_material = material.clone()
        original_max_z = get_atomic_coordinates_extremum(new_material, use_cartesian_coordinates=True)
        material_with_additional_layers = self.create_material_with_additional_layers(new_material, thickness)
        added_layers_max_z = get_atomic_coordinates_extremum(
            material_with_additional_layers, use_cartesian_coordinates=True
        )
        if condition is None:
            condition = self._default_condition

        atoms_within_island = filter_by_condition_on_coordinates(
            material=material_with_additional_layers,
            condition=condition,
            use_cartesian_coordinates=use_cartesian_coordinates,
        )
        # Filter atoms in the added layers between the original and added layers
        island_material = filter_by_box(
            material=atoms_within_island,
            min_coordinate=[0, 0, original_max_z],
            max_coordinate=[material.lattice.a, material.lattice.b, added_layers_max_z],
            use_cartesian_coordinates=True,
        )

        return self.merge_slab_and_defect(island_material, new_material)

    def _generate(self, configuration: _ConfigurationType) -> List[_GeneratedItemType]:
        condition_callable = configuration.condition.condition
        return [
            self.create_island(
                material=configuration.crystal,
                condition=condition_callable,
                thickness=configuration.number_of_added_layers,
                use_cartesian_coordinates=configuration.use_cartesian_coordinates,
            )
        ]
