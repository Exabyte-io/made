from typing import List, Tuple, Union

import numpy as np
from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from pydantic import BaseModel, Field

from mat3ra.made.material import Material

from .basis import BasisMaterialAnalyzer
from ..build_components.metadata.material_with_build_metadata import MaterialWithBuildMetadata
from ..modify import wrap_to_unit_cell


class LatticeSwapParameters(BaseModel):
    model_config = {"arbitrary_types_allowed": True}

    permutation: List[Tuple[int, int]] = Field(..., description="List of (source_idx, sign) tuples for each new vector")


class LatticeSwapDetectionResult(LatticeSwapParameters):
    is_swapped: bool = Field(..., description="Whether a lattice swap was detected")
    permutation: Matrix3x3Schema = Field(..., description="Transformation matrix representing the swap")
    new_lattice_vectors: List = Field(..., description="New lattice vectors")
    new_coordinates: List = Field(..., description="New lattice coordinates")
    confidence: float = Field(..., description="Confidence score (0.0 to 1.0)")


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

        possible_transformation_matrices = [
            # np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]]) / 2,
            # np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]]) / 2,
            # np.array([[1, 1, 0], [-1, 1, 0], [0, 0, 2]]) / 2,
            # np.array([[1, -1, 0], [1, 1, 0], [0, 0, 2]]) / 2,
            # Direct swaps
            np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]]),
            np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]]),
            np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]]),
            # Inverted swaps
            np.array([[0, -1, 0], [-1, 0, 0], [0, 0, 1]]),
            np.array([[0, 0, -1], [0, 1, 0], [-1, 0, 0]]),
            np.array([[1, 0, 0], [0, 0, -1], [0, -1, 0]]),
        ]

        max_score = 0
        best_match = None

        lattice_vectors = original_material.lattice.vector_arrays
        coordinates = original_material.basis.coordinates.values
        original_analyzer = BasisMaterialAnalyzer(material=original_material)
        original_fingerprint = original_analyzer.get_layer_fingerprint(layer_thickness)

        for matrix in possible_transformation_matrices:
            new_lattice_vectors = matrix @ lattice_vectors
            new_coordinates = np.linalg.inv(matrix) @ np.array(coordinates).T
            new_coordinates = new_coordinates.T.tolist()

            transformed_material = original_material.clone()
            transformed_material.set_lattice_vectors_from_array(new_lattice_vectors.tolist())
            transformed_material.basis.coordinates.values = new_coordinates
            transformed_material = wrap_to_unit_cell(transformed_material)

            new_analyzer = BasisMaterialAnalyzer(material=transformed_material)
            new_fingerprint = new_analyzer.get_layer_fingerprint(layer_thickness)
            score = original_fingerprint.get_similarity_score(new_fingerprint)

            if score > max_score:
                best_match = LatticeSwapDetectionResult(
                    is_swapped=True,
                    permutation=Matrix3x3Schema(root=matrix.tolist()),
                    new_lattice_vectors=new_lattice_vectors,
                    new_coordinates=new_coordinates,
                    confidence=score,
                )
                max_score = score

        if best_match and best_match.confidence >= threshold:
            return best_match

        return LatticeSwapDetectionResult(
            is_swapped=False,
            permutation=Matrix3x3Schema(root=np.eye(3).tolist()),
            new_lattice_vectors=self.material.lattice.vector_arrays,
            new_coordinates=self.material.basis.coordinates.values,
            confidence=0.0,
        )

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
            corrected_material = target.clone()
            corrected_material.set_lattice_vectors_from_array(swap_info.new_lattice_vectors)
            corrected_material.basis.coordinates.values = swap_info.new_coordinates
            return corrected_material
        return self.material
