from typing import List, Callable, Optional, Union

import numpy as np

from mat3ra.made.material import Material
from .configuration import (
    PointDefectConfigurationLegacy,
    AdatomSlabPointDefectConfiguration,
    IslandSlabDefectConfiguration,
    TerraceSlabDefectConfiguration,
    PointDefectPairConfiguration,
)
from .factories import DefectBuilderFactory
from .point.builders import PointDefectBuilder
from .slab.builders import SlabDefectBuilder
from ...analyze.other import (
    get_atomic_coordinates_extremum,
)
from ...build import BaseBuilder
from ...modify import (
    filter_by_box,
    filter_by_condition_on_coordinates,
    translate_to_z_level,
    rotate,
)
from ...utils import (
    coordinate as CoordinateCondition,
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


#
#
# class IslandSlabDefectBuilder(SlabDefectBuilder):
#     """
#     Legacy builder for island defects. Now uses the new pattern with VoidSite and IslandDefectConfiguration.
#     """
#     _ConfigurationType: type(IslandSlabDefectConfiguration) = IslandSlabDefectConfiguration  # type: ignore
#     _GeneratedItemType: Material = Material
#
#     def _generate(self, configuration: _ConfigurationType) -> List[_GeneratedItemType]:
#         from mat3ra.made.tools.build.defect.slab.helpers import create_island_defect
#
#         return [
#             create_island_defect(
#                 slab=configuration.crystal,
#                 condition=configuration.condition,
#                 number_of_added_layers=configuration.number_of_added_layers,
#             )
#         ]


class TerraceSlabDefectBuilder(SlabDefectBuilder):
    _ConfigurationType: type(TerraceSlabDefectConfiguration) = TerraceSlabDefectConfiguration  # type: ignore
    _GeneratedItemType: Material = Material

    def _calculate_cut_direction_vector(self, material: Material, cut_direction: List[int]):
        """
        Calculate the normalized cut direction vector in Cartesian coordinates.

        Args:
            material: The material to get the lattice vectors from.
            cut_direction: The direction of the cut in lattice directions.

        Returns:
            The normalized cut direction vector in Cartesian coordinates.
        """
        np_cut_direction = np.array(cut_direction)
        direction_vector = np.dot(np.array(material.basis.cell.vector_arrays), np_cut_direction)
        normalized_direction_vector = direction_vector / np.linalg.norm(direction_vector)
        return normalized_direction_vector

    def _calculate_height_cartesian(self, material: Material, new_material: Material):
        """
        Calculate the height of the added layers in Cartesian coordinates.

        Args:
            material: The original material.
            new_material: The material with the added layers.

        Returns:
            The height of the added layers in Cartesian coordinates.
        """
        original_max_z = get_atomic_coordinates_extremum(material, use_cartesian_coordinates=True)
        added_layers_max_z = get_atomic_coordinates_extremum(new_material, use_cartesian_coordinates=True)
        height_cartesian = added_layers_max_z - original_max_z
        return height_cartesian

    def _calculate_rotation_parameters(
        self, original_material: Material, new_material: Material, normalized_direction_vector: List[float]
    ):
        """
        Calculate the necessary rotation angle and axis.

        Args:
            original_material: The original material.
            new_material: The material with the added layers.
            normalized_direction_vector: The normalized cut direction vector in Cartesian coordinates.

        Returns:
            The rotation angle, normalized rotation axis, and delta length.
        """
        height_cartesian = self._calculate_height_cartesian(original_material, new_material)
        cut_direction_xy_proj_cart = np.linalg.norm(
            np.dot(np.array(new_material.lattice.vector_arrays), normalized_direction_vector)
        )
        # Slope of the terrace along the cut direction
        hypotenuse = np.linalg.norm([height_cartesian, cut_direction_xy_proj_cart])
        angle = np.arctan(height_cartesian / cut_direction_xy_proj_cart) * 180 / np.pi
        normalized_rotation_axis = np.cross(normalized_direction_vector, [0, 0, 1]).tolist()
        delta_length = hypotenuse - cut_direction_xy_proj_cart
        return angle, normalized_rotation_axis, delta_length

    def _increase_lattice_size(
        self, material: Material, length_increase: float, direction_of_increase: List[float]
    ) -> Material:
        """
        Increase the lattice size in a specific direction.

        When the material is rotated to maintain periodic boundary conditions (PBC),
        the periodicity in the X and Y directions changes.
        Therefore, the lattice size must be increased to fit the new structure dimensions.

        If the terrace plane is normal to the Z direction, it becomes larger than the previous XY plane of the material
        because it forms a hypotenuse between PBC points.
        This method adjusts the lattice vectors to accommodate this change.

        Args:
            material: The material to increase the lattice size of.
            length_increase: The increase in length.
            direction_of_increase: The direction of the increase.

        Returns:
            The material with the increased lattice size.
        """
        vector_a, vector_b = material.lattice.vectors.a, material.lattice.vectors.b
        norm_a, norm_b = vector_a.norm, vector_b.norm

        delta_a_cart = np.dot(vector_a.value, np.array(direction_of_increase)) * length_increase / norm_a
        delta_b_cart = np.dot(vector_b.value, np.array(direction_of_increase)) * length_increase / norm_b

        scaling_matrix = np.eye(3)
        scaling_matrix[0, 0] += delta_a_cart / norm_a
        scaling_matrix[1, 1] += delta_b_cart / norm_b

        new_lattice = material.lattice.get_scaled_by_matrix(scaling_matrix.tolist())
        material.set_lattice(new_lattice)
        return material

    def _update_material_name(self, material: Material, configuration: _ConfigurationType) -> Material:
        new_material = super()._update_material_name(material, configuration)
        new_name = (
            f"{new_material.name}, {configuration.number_of_added_layers}-step Terrace {configuration.cut_direction}"
        )
        new_material.name = new_name
        return new_material

    def create_terrace(
        self,
        material: Material,
        cut_direction: Optional[List[int]] = None,
        pivot_coordinate: Optional[List[float]] = None,
        number_of_added_layers: int = 1,
        use_cartesian_coordinates: bool = False,
        rotate_to_match_pbc: bool = True,
    ) -> Material:
        """
        Create a terrace at the specified position on the surface of the material.

        Args:
            material: The material to add the terrace to.
            cut_direction: The direction of the cut in lattice directions.
            pivot_coordinate: The center position of the terrace.
            number_of_added_layers: The number of added layers to the slab which will form the terrace
            use_cartesian_coordinates: Whether to use Cartesian coordinates for the center position.
            rotate_to_match_pbc: Whether to rotate the material to match the periodic boundary conditions.
        Returns:
            The material with the terrace added.
        """
        if cut_direction is None:
            cut_direction = [0, 0, 1]
        if pivot_coordinate is None:
            pivot_coordinate = [0.5, 0.5, 0.5]

        new_material = material.clone()
        material_with_additional_layers = self.create_material_with_additional_layers(
            new_material, number_of_added_layers
        )

        normalized_direction_vector = self._calculate_cut_direction_vector(material, cut_direction)
        condition = CoordinateCondition.PlaneCoordinateCondition(
            plane_normal=normalized_direction_vector,
            plane_point_coordinate=pivot_coordinate,
        ).condition
        atoms_within_terrace = filter_by_condition_on_coordinates(
            material=material_with_additional_layers,
            condition=condition,
            use_cartesian_coordinates=use_cartesian_coordinates,
        )
        merged_material = self.merge_slab_and_defect(new_material, atoms_within_terrace)

        angle, normalized_rotation_axis, delta_length = self._calculate_rotation_parameters(
            material, merged_material, normalized_direction_vector
        )
        result_material = translate_to_z_level(merged_material, "center")

        if rotate_to_match_pbc:
            adjusted_material = self._increase_lattice_size(result_material, delta_length, normalized_direction_vector)
            result_material = rotate(material=adjusted_material, axis=normalized_rotation_axis, angle=angle)
        return result_material

    def _generate(self, configuration: _ConfigurationType) -> List[_GeneratedItemType]:
        return [
            self.create_terrace(
                material=configuration.crystal,
                cut_direction=configuration.cut_direction,
                pivot_coordinate=configuration.pivot_coordinate,
                number_of_added_layers=configuration.number_of_added_layers,
                use_cartesian_coordinates=configuration.use_cartesian_coordinates,
            )
        ]
