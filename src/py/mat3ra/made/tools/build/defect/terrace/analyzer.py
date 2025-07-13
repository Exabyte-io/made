from typing import List

import numpy as np

from mat3ra.made.tools.analyze.slab import SlabMaterialAnalyzer


class TerraceMaterialAnalyzer(SlabMaterialAnalyzer):
    def calculate_cut_direction_vector(self, cut_direction: List[int]):
        np_cut_direction = np.array(cut_direction)
        direction_vector = np.dot(np.array(self.material.basis.cell.vector_arrays), np_cut_direction)
        normalized_direction_vector = direction_vector / np.linalg.norm(direction_vector)
        return normalized_direction_vector

    def get_new_lattice_vectors(self, length_increase: float, direction_of_increase: List[float]) -> List[float]:
        """
        Increase the lattice size in a specific direction.

        When the material is rotated to maintain periodic boundary conditions (PBC),
        the periodicity in the X and Y directions changes.
        Therefore, the lattice size must be increased to fit the new structure dimensions.

        If the terrace plane is normal to the Z direction, it becomes larger than the previous XY plane of the material
        because it forms a hypotenuse between PBC points.
        This method adjusts the lattice vectors to accommodate this change.

        Args:
            length_increase: The increase in length.
            direction_of_increase: The direction of the increase.

        Returns:
            The material with the increased lattice size.
        """
        vector_a, vector_b = self.material.lattice.vectors.a, self.material.lattice.vectors.b
        norm_a, norm_b = vector_a.norm, vector_b.norm

        delta_a_cart = np.dot(vector_a.value, np.array(direction_of_increase)) * length_increase / norm_a
        delta_b_cart = np.dot(vector_b.value, np.array(direction_of_increase)) * length_increase / norm_b

        scaling_matrix = np.eye(3)
        scaling_matrix[0, 0] += delta_a_cart / norm_a
        scaling_matrix[1, 1] += delta_b_cart / norm_b

        new_material = self.material.clone()
        new_lattice = new_material.lattice.get_scaled_by_matrix(scaling_matrix.tolist())
        return new_lattice.vector_arrays.tolist()

    # def calculate_height_cartesian(self, new_material: MaterialWithBuildMetadata):
    #     """
    #     Calculate the height of the added layers in Cartesian coordinates.
    #
    #     Args:
    #         material: The original material.
    #         new_material: The material with the added layers.
    #
    #     Returns:
    #         The height of the added layers in Cartesian coordinates.
    #     """
    #     original_max_z = get_atomic_coordinates_extremum(material, use_cartesian_coordinates=True)
    #     added_layers_max_z = get_atomic_coordinates_extremum(new_material, use_cartesian_coordinates=True)
    #     height_cartesian = added_layers_max_z - original_max_z
    #     return height_cartesian
    #
    # def _calculate_rotation_parameters(
    #     self,
    #     original_material: MaterialWithBuildMetadata,
    #     new_material: MaterialWithBuildMetadata,
    #     normalized_direction_vector: List[float],
    # ):
    #     """
    #     Calculate the necessary rotation angle and axis.
    #
    #     Args:
    #         original_material: The original material.
    #         new_material: The material with the added layers.
    #         normalized_direction_vector: The normalized cut direction vector in Cartesian coordinates.
    #
    #     Returns:
    #         The rotation angle, normalized rotation axis, and delta length.
    #     """
    #     height_cartesian = self._calculate_height_cartesian(original_material, new_material)
    #     cut_direction_xy_proj_cart = np.linalg.norm(
    #         np.dot(np.array(new_material.lattice.vector_arrays), normalized_direction_vector)
    #     )
    #     # Slope of the terrace along the cut direction
    #     hypotenuse = np.linalg.norm([height_cartesian, cut_direction_xy_proj_cart])
    #     angle = np.arctan(height_cartesian / cut_direction_xy_proj_cart) * 180 / np.pi
    #     normalized_rotation_axis = np.cross(normalized_direction_vector, [0, 0, 1]).tolist()
    #     delta_length = hypotenuse - cut_direction_xy_proj_cart
    #     return angle, normalized_rotation_axis, delta_length
    #
