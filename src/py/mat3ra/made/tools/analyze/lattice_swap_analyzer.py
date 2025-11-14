from typing import List, Union

import numpy as np
from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from pydantic import BaseModel, Field

from mat3ra.made.material import Material
from .basis import BasisMaterialAnalyzer
from ..build_components.metadata.material_with_build_metadata import MaterialWithBuildMetadata
from ..modify import wrap_to_unit_cell
from ...lattice import Lattice

POSSIBLE_TRANSFORMATION_MATRICES = [
    [[1, 0, 0], [0, 1, 0], [0, 0, 1]],  # Identity - no transformation
    # Direct swaps: new[i] = old[j]
    [[0, 1, 0], [1, 0, 0], [0, 0, 1]],
    [[0, 0, 1], [0, 1, 0], [1, 0, 0]],
    [[1, 0, 0], [0, 0, 1], [0, 1, 0]],
    # Swaps with sign flips
    [[0, 0, 1], [1, 0, 0], [0, -1, 0]],
    [[0, 0, -1], [1, 0, 0], [0, 1, 0]],
    [[0, 1, 0], [0, 0, 1], [-1, 0, 0]],
    [[0, -1, 0], [0, 0, 1], [1, 0, 0]],
    # Inverted swaps
    [[0, -1, 0], [-1, 0, 0], [0, 0, 1]],
    [[0, 0, -1], [0, 1, 0], [-1, 0, 0]],
    [[1, 0, 0], [0, 0, -1], [0, -1, 0]],
    # Rotations around x-axis: -90 and +90 degrees
    [[1, 0, 0], [0, 0, 1], [0, -1, 0]],  # -90° around x
    [[1, 0, 0], [0, 0, -1], [0, 1, 0]],  # +90° around x
    # Rotations around y-axis: -90 and +90 degrees
    [[0, 0, -1], [0, 1, 0], [1, 0, 0]],  # -90° around y
    [[0, 0, 1], [0, 1, 0], [-1, 0, 0]],  # +90° around y
    # Rotations around z-axis: -90 and +90 degrees
    [[0, 1, 0], [-1, 0, 0], [0, 0, 1]],  # -90° around z
    [[0, -1, 0], [1, 0, 0], [0, 0, 1]],  # +90° around z
    # 180° rotations around axes
    [[1, 0, 0], [0, -1, 0], [0, 0, -1]],  # 180° around x
    [[-1, 0, 0], [0, 1, 0], [0, 0, -1]],  # 180° around y
    [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],  # 180° around z
    # Mirrors (reflections)
    [[1, 0, 0], [0, 1, 0], [0, 0, -1]],  # Mirror along z
    [[1, 0, 0], [0, -1, 0], [0, 0, 1]],  # Mirror along y
    [[-1, 0, 0], [0, 1, 0], [0, 0, 1]],  # Mirror along x
]


class LatticeSwapDetectionResult(BaseModel):
    model_config = {"arbitrary_types_allowed": True}
    is_swapped: bool = Field(..., description="Whether a lattice swap was detected")
    permutation: Matrix3x3Schema = Field(..., description="Transformation matrix representing the swap")
    new_lattice_vectors: List = Field(..., description="New lattice vectors")
    new_coordinates: List = Field(..., description="New lattice coordinates")
    confidence: float = Field(..., description="Confidence score (0.0 to 1.0)")

    @property
    def new_lattice(self) -> Lattice:
        return Lattice.from_vectors_array(vectors=self.new_lattice_vectors)


class MaterialLatticeSwapAnalyzer(BaseModel):
    """
    Analyzer to detect lattice vector swaps/permutations between materials.

    This detects when lattice vectors have been reoriented (e.g., a->a, b->c, c->-b)
    rather than the basis being rotated within the lattice.
    """

    material: Union[Material, MaterialWithBuildMetadata]
    tolerance: float = 0.01

    def detect_swap_from_original(
        self, original_material: MaterialWithBuildMetadata, layer_thickness: float = 1.0, threshold: float = 0.1
    ) -> LatticeSwapDetectionResult:
        """
        Detect lattice swap from the original material.

        We apply every transformation matrix to lattice and basis, compute the fingerprint along z, then compare
        to the original material's fingerprint. The best matching transformation (if above threshold) is considered
        a detected swap.


        Args:
            original_material: The original material before transformation
            layer_thickness: Thickness of layers for fingerprint comparison
            threshold: Minimum improvement threshold to consider a swap detected

        Returns:
            LatticeSwapDetectionResult: Swap detection results with permutation and new lattice
        """
        possible_transformation_matrices = list(map(np.array, POSSIBLE_TRANSFORMATION_MATRICES))

        max_score = 0
        best_match = None

        target_analyzer = BasisMaterialAnalyzer(material=original_material)
        target_fingerprint = target_analyzer.get_layer_fingerprint(layer_thickness)

        current_lattice_vectors = np.array(self.material.lattice.vector_arrays)
        current_coordinates = self.material.basis.coordinates.values

        for matrix in possible_transformation_matrices:
            new_lattice_vectors = (matrix @ current_lattice_vectors.T).tolist()
            new_coordinates = (np.linalg.inv(matrix) @ np.array(current_coordinates).T).T.tolist()

            transformed_material = self.material.clone()
            transformed_material.set_lattice_vectors_from_array(new_lattice_vectors)
            transformed_material.basis.coordinates.values = new_coordinates
            transformed_material = wrap_to_unit_cell(transformed_material)

            new_analyzer = BasisMaterialAnalyzer(material=transformed_material)
            new_fingerprint = new_analyzer.get_layer_fingerprint(layer_thickness)
            score = target_fingerprint.get_similarity_score(new_fingerprint)

            is_identity = np.allclose(matrix, np.eye(3), atol=self.tolerance)

            if score > max_score:
                best_match = LatticeSwapDetectionResult(
                    is_swapped=not is_identity,
                    permutation=Matrix3x3Schema(root=matrix.tolist()),
                    new_lattice_vectors=(
                        self.material.lattice.vector_arrays if not is_identity else new_lattice_vectors
                    ),
                    new_coordinates=(self.material.basis.coordinates.values if not is_identity else new_coordinates),
                    confidence=score,
                )
                max_score = score

        if best_match and best_match.is_swapped and best_match.confidence >= threshold:
            return best_match

        return best_match

    def correct_material_to_match_target(
        self,
        target: MaterialWithBuildMetadata,
        layer_thickness: float = 1.0,
        threshold: float = 0.1,
    ) -> MaterialWithBuildMetadata:
        """
        Correct the material's lattice to match the target material's orientation if a swap is detected.

        Args:
            target: The target material to match
            layer_thickness: Thickness of layers for fingerprint comparison
            threshold: Minimum improvement threshold to consider a swap detected

        Returns:
            MaterialWithBuildMetadata: Corrected material with lattice matching the target orientation
        """
        swap_info = self.detect_swap_from_original(target, layer_thickness, threshold)
        if swap_info.is_swapped:
            # Apply the detected transformation to self.material to match target
            matrix_list = [list(row.root) if hasattr(row, "root") else row for row in swap_info.permutation.root]
            matrix_array = np.array(matrix_list)
            current_lattice_vectors = np.array(self.material.lattice.vector_arrays)
            current_coordinates = self.material.basis.coordinates.values

            corrected_lattice_vectors = (matrix_array @ current_lattice_vectors.T).tolist()
            corrected_coordinates = (np.linalg.inv(matrix_array) @ np.array(current_coordinates).T).T.tolist()

            corrected_material = self.material.clone()
            corrected_material.set_lattice_vectors_from_array(corrected_lattice_vectors)
            corrected_material.basis.coordinates.values = corrected_coordinates
            corrected_material = wrap_to_unit_cell(corrected_material)
            return corrected_material
        return self.material
