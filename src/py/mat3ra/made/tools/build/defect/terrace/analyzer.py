from typing import List
import numpy as np

from mat3ra.made.tools.analyze.slab import SlabMaterialAnalyzer


class TerraceMaterialAnalyzer(SlabMaterialAnalyzer):
    cut_direction: List[int] = [0, 0, 1]
    number_of_added_layers: float = 1.0

    @property
    def _lattice_matrix(self) -> np.ndarray:
        return np.array(self.material.lattice.vector_arrays)

    @property
    def cut_direction_vector(self) -> np.ndarray:
        direction = np.dot(self._lattice_matrix, np.array(self.cut_direction))
        return direction / np.linalg.norm(direction)

    @property
    def _height_cartesian(self) -> float:
        return self.layer_thickness * self.number_of_added_layers

    @property
    def _cut_direction_xy_proj_cart(self) -> float:
        return np.linalg.norm(np.dot(self._lattice_matrix, self.cut_direction_vector))

    @property
    def rotation_axis(self) -> List[int]:
        axis = np.cross(np.array(self.cut_direction), [0, 0, 1])
        return np.round(axis / np.linalg.norm(axis)).astype(int).tolist()

    @property
    def angle(self) -> float:
        return np.degrees(np.arctan(self._height_cartesian / self._cut_direction_xy_proj_cart))

    @property
    def _length_increase(self) -> float:
        return (
            np.linalg.norm([self._height_cartesian, self._cut_direction_xy_proj_cart])
            - self._cut_direction_xy_proj_cart
        )

    @property
    def new_lattice_vectors(self) -> np.ndarray:
        dir_vec = self.cut_direction_vector
        length_inc = self._length_increase

        vector_a = self.material.lattice.vectors.a
        vector_b = self.material.lattice.vectors.b

        norm_a, norm_b = vector_a.norm, vector_b.norm
        delta_a_cart = np.dot(vector_a.value, dir_vec) * length_inc / norm_a
        delta_b_cart = np.dot(vector_b.value, dir_vec) * length_inc / norm_b

        scaling_matrix = np.eye(3)
        scaling_matrix[0, 0] += delta_a_cart / norm_a
        scaling_matrix[1, 1] += delta_b_cart / norm_b

        new_material = self.material.clone()
        new_lattice = new_material.lattice.get_scaled_by_matrix(scaling_matrix.tolist())
        return new_lattice.vector_arrays
