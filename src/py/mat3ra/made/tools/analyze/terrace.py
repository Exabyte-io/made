from typing import List

import numpy as np

from .slab import SlabMaterialAnalyzer


class TerraceMaterialAnalyzer(SlabMaterialAnalyzer):
    """
    Used to calculate the adjustment values to make a stepped-surface (terrace) slab
    periodic.

    Build sequence
    --------------
    1. **Add extra layers** on top of the parent slab.
    2. **Cut** those layers with the crystallographic plane given by
       ``cut_direction`` → isolates the terrace step.
    3. **Stack** the step onto the original slab → creates a terrace.
    4. **Stretch lattice** in-plane (``new_lattice_vectors``) so the rotation can be performed inside the unit cell.
    5. **Rotate** the whole slab by ``angle`` around ``rotation_axis`` so the
       crystal structure is matched with the repetition in the PBC.

    Key outputs
    -----------
    cut_direction_vector : unit Cartesian vector normal to the cut plane
    rotation_axis        : axis to rotate the slab around
    angle                : rotation angle in degrees
    new_lattice_vectors  : stretched lattice vectors used before rotation
    """

    cut_direction: List[int] = [0, 0, 1]
    number_of_added_layers: float = 1.0

    @property
    def _lattice_matrix(self) -> np.ndarray:
        return np.array(self.material.lattice.vector_arrays)

    @property
    def _height_cartesian(self) -> float:
        return self.layer_thickness * self.number_of_added_layers

    @property
    def _cut_direction_xy_proj_cart(self) -> float:
        return np.linalg.norm(np.dot(self._lattice_matrix, self.cut_direction_vector))

    @property
    def _length_increment(self) -> float:
        return (
            np.linalg.norm([self._height_cartesian, self._cut_direction_xy_proj_cart])
            - self._cut_direction_xy_proj_cart
        )

    @property
    def cut_direction_vector(self) -> np.ndarray:
        """
        Get the normalized Cartesian vector normal to the terrace cut plane.

        The terrace is constructed by cutting the material with a crystallographic
        plane defined by `cut_direction`
        """
        direction = np.dot(self._lattice_matrix, np.array(self.cut_direction))
        return direction / np.linalg.norm(direction)

    @property
    def rotation_axis(self) -> List[int]:
        """
        Calculate the axis of rotation to align the terrace with its repetition in the PBC.
        """
        axis = np.cross(np.array(self.cut_direction), [0, 0, 1])
        return np.round(axis / np.linalg.norm(axis)).astype(int).tolist()

    @property
    def angle(self) -> float:
        """
        Calculate the angle of the tilt needed to align the terrace with its repetition in the PBC.
        """
        return np.degrees(np.arctan(self._height_cartesian / self._cut_direction_xy_proj_cart))

    @property
    def new_lattice_vectors(self) -> List[List[float]]:
        """
        Compute updated lattice vectors to accommodate terrace tilt before rotation.

        The terrace step introduces an additional height, so to match it across PBC,
        the in-plane lattice vectors must be stretched. This method calculates the
        additional in-plane length needed and applies directional scaling to both
        a and b vectors, in proportion to their projection along the terrace normal.
        """
        direction_vector = self.cut_direction_vector
        length_increment = self._length_increment

        vector_a = self.material.lattice.vectors.a
        vector_b = self.material.lattice.vectors.b

        norm_a, norm_b = vector_a.norm, vector_b.norm
        delta_a_cart = np.dot(vector_a.value, direction_vector) * length_increment / norm_a
        delta_b_cart = np.dot(vector_b.value, direction_vector) * length_increment / norm_b

        scaling_matrix = np.eye(3)
        scaling_matrix[0, 0] += delta_a_cart / norm_a
        scaling_matrix[1, 1] += delta_b_cart / norm_b

        new_material = self.material.clone()
        new_lattice = new_material.lattice.get_scaled_by_matrix(scaling_matrix.tolist())
        return new_lattice.vector_arrays
