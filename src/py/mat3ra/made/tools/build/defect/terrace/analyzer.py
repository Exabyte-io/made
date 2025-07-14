from typing import List

import numpy as np

from mat3ra.made.tools.analyze.slab import SlabMaterialAnalyzer
from mat3ra.made.tools.build import MaterialWithBuildMetadata
from mat3ra.made.utils import get_atomic_coordinates_extremum


class TerraceMaterialAnalyzer(SlabMaterialAnalyzer):
    def calculate_cut_direction_vector(self, cut_direction: List[int]):
        np_cut_direction = np.array(cut_direction)
        direction_vector = np.dot(np.array(self.material.basis.cell.vector_arrays), np_cut_direction)
        normalized_direction_vector = direction_vector / np.linalg.norm(direction_vector)
        return normalized_direction_vector

    def get_rotation_axis(self, cut_direction: List[int]) -> List[int]:
        normalized_rotation_axis = np.cross(cut_direction, [0, 0, 1])
        norm = np.linalg.norm(normalized_rotation_axis)
        normalized = normalized_rotation_axis / norm
        return np.round(normalized).astype(int).tolist()

    def get_angle(self, cut_direction, number_of_added_layers: float) -> float:
        """
        Calculate the angle of rotation based on the number of added layers.

        Args:
            cut_direction: The direction of the cut in Miller indices.
            number_of_added_layers: The number of layers added to the slab.

        Returns:
            The angle in degrees.
        """
        height_cartesian = self.layer_thickness * number_of_added_layers
        cut_direction_xy_proj_cart = np.linalg.norm(
            np.dot(np.array(self.material.lattice.vector_arrays), self.calculate_cut_direction_vector(cut_direction))
        )
        angle = np.arctan(height_cartesian / cut_direction_xy_proj_cart) * 180 / np.pi
        return angle

    def get_length_increase(self, number_of_added_layers: float) -> float:
        """
        Calculate the increase in length based on the number of added layers.

        Args:
            number_of_added_layers: The number of layers added to the slab.

        Returns:
            The increase in length.
        """
        height_cartesian = self.layer_thickness * number_of_added_layers
        cut_direction_xy_proj_cart = np.linalg.norm(
            np.dot(np.array(self.material.lattice.vector_arrays), self.calculate_cut_direction_vector([0, 0, 1]))
        )
        hypotenuse = np.linalg.norm([height_cartesian, cut_direction_xy_proj_cart])
        delta_length = hypotenuse - cut_direction_xy_proj_cart
        return delta_length

    def get_new_lattice_vectors(self, cut_direction, number_of_added_layers: float) -> List[List[float]]:
        """
        Increase the lattice size in a specific direction.

        When the material is rotated to maintain periodic boundary conditions (PBC),
        the periodicity in the X and Y directions changes.
        Therefore, the lattice size must be increased to fit the new structure dimensions.

        If the terrace plane is normal to the Z direction, it becomes larger than the previous XY plane of the material
        because it forms a hypotenuse between PBC points.
        This method adjusts the lattice vectors to accommodate this change.

        Args:
            cut_direction: The direction of the cut in Miller indices.
            number_of_added_layers: The number of layers added to the slab.

        Returns:
            The material with the increased lattice size.
        """

        direction_of_increase = self.calculate_cut_direction_vector(cut_direction)
        length_increase = self.get_length_increase(number_of_added_layers)
        vector_a, vector_b = self.material.lattice.vectors.a, self.material.lattice.vectors.b
        norm_a, norm_b = vector_a.norm, vector_b.norm

        delta_a_cart = np.dot(vector_a.value, np.array(direction_of_increase)) * length_increase / norm_a
        delta_b_cart = np.dot(vector_b.value, np.array(direction_of_increase)) * length_increase / norm_b

        scaling_matrix = np.eye(3)
        scaling_matrix[0, 0] += delta_a_cart / norm_a
        scaling_matrix[1, 1] += delta_b_cart / norm_b

        new_material = self.material.clone()
        new_lattice = new_material.lattice.get_scaled_by_matrix(scaling_matrix.tolist())
        return new_lattice.vector_arrays

    def calculate_height_cartesian(self, new_material: MaterialWithBuildMetadata):
        """
        Calculate the height of the added layers in Cartesian coordinates.

        Args:
            new_material: The material with the added layers.

        Returns:
            The height of the added layers in Cartesian coordinates.
        """
        original_max_z = get_atomic_coordinates_extremum(self.material, use_cartesian_coordinates=True)
        added_layers_max_z = get_atomic_coordinates_extremum(new_material, use_cartesian_coordinates=True)
        height_cartesian = added_layers_max_z - original_max_z
        return height_cartesian
