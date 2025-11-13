from typing import List, Tuple, Union

import numpy as np
from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.rotation_analyzer import MaterialRotationAnalyzer
from mat3ra.made.tools.modify import wrap_to_unit_cell, translate_to_z_level
from pydantic import BaseModel, Field

from ..build_components.metadata.material_with_build_metadata import MaterialWithBuildMetadata
from ..operations.core.unary import rotate
from ...lattice import Lattice


class LatticeSwapParameters(BaseModel):
    model_config = {"arbitrary_types_allowed": True}

    permutation: List[Tuple[int, int]] = Field(..., description="List of (source_idx, sign) tuples for each new vector")


class LatticeSwapDetectionResult(LatticeSwapParameters):
    is_swapped: bool = Field(..., description="Whether a lattice swap was detected")
    permutation: List[Tuple[int, int]] = Field(..., description="Permutation that was detected")
    new_lattice: Lattice = Field(..., description="New lattice configuration with swapped parameters")
    confidence: float = Field(..., description="Confidence score (0.0 to 1.0)")


class MaterialLatticeSwapAnalyzer(MaterialRotationAnalyzer):
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
        Detect lattice vector swap/permutation of the current material compared to the original material.

        Args:
            original_material: The original material before transformation
            layer_thickness: Thickness of layers for fingerprint comparison
            threshold: Minimum improvement threshold to consider a swap detected

        Returns:
            LatticeSwapDetectionResult: Swap detection results with permutation and new lattice
        """

        original_lattice = original_material.lattice
        current_lattice = self.material.lattice

        rotation_result = self.detect_rotation_from_original(original_material, layer_thickness, threshold)

        if not rotation_result.is_rotated:
            return LatticeSwapDetectionResult(
                is_swapped=False,
                permutation=[(0, 1), (1, 1), (2, 1)],
                new_lattice=original_lattice,
                confidence=rotation_result.confidence,
            )

        original_vectors = np.array(original_lattice.vector_arrays)
        current_vectors = np.array(current_lattice.vector_arrays)

        if np.allclose(original_vectors, current_vectors, atol=self.tolerance):
            return LatticeSwapDetectionResult(
                is_swapped=True,
                permutation=[(0, 1), (1, 1), (2, 1)],
                new_lattice=current_lattice,
                confidence=rotation_result.confidence,
            )

        permutation = self._rotation_matrix_to_permutation(rotation_result.rotation_matrix)
        new_lattice = self._apply_permutation_to_lattice(original_material.lattice, permutation)

        return LatticeSwapDetectionResult(
            is_swapped=True,
            permutation=permutation,
            new_lattice=new_lattice,
            confidence=rotation_result.confidence,
        )

    def _rotation_matrix_to_permutation(self, rotation_matrix: np.ndarray) -> List[Tuple[int, int]]:
        """
        Convert a rotation matrix to a lattice vector permutation.

        The rotation matrix R transforms coordinates. For lattice vectors, we need to determine
        which original vector (0=a, 1=b, 2=c) maps to each new position, and whether there's a sign flip.

        Args:
            rotation_matrix: 3x3 rotation matrix

        Returns:
            List of (source_idx, sign) tuples where source_idx is the original vector index (0,1,2)
            and sign is 1 or -1 indicating whether the vector should be flipped
        """
        permutation = []
        used_indices = set()

        for new_idx in range(3):
            col = rotation_matrix[:, new_idx]
            best_match_idx = None
            best_match_value = 0.0

            for orig_idx in range(3):
                if orig_idx in used_indices:
                    continue
                value = abs(col[orig_idx])
                if value > best_match_value:
                    best_match_value = value
                    best_match_idx = orig_idx

            if best_match_idx is not None and best_match_value > self.tolerance:
                sign = 1 if col[best_match_idx] >= 0 else -1
                permutation.append((best_match_idx, sign))
                used_indices.add(best_match_idx)
            else:
                permutation.append((new_idx, 1))

        return permutation

    def _apply_permutation_to_lattice(self, original_lattice: Lattice, permutation: List[Tuple[int, int]]) -> Lattice:
        vectors = original_lattice.vector_arrays
        new_vectors = []

        for source_idx, sign in permutation:
            vector = np.array(vectors[source_idx]) * sign
            new_vectors.append(vector.tolist())

        new_lattice = Lattice.from_vectors_array(
            vectors=new_vectors, units=original_lattice.units, type=original_lattice.type
        )
        return new_lattice

    def _invert_permutation(self, permutation: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
        inverse = [None] * 3
        for new_idx, (source_idx, sign) in enumerate(permutation):
            inverse[source_idx] = (new_idx, sign)
        return inverse

    def correct_lattice_to_match_original(
        self, original_material: MaterialWithBuildMetadata, layer_thickness: float = 1.0, threshold: float = 0.1
    ) -> Lattice:
        """
        Correct the current material's lattice to match the original material's orientation.

        This method detects any rotation/swap between the current and original material,
        and applies the inverse transformation to correct the current material's lattice
        back to the original orientation.

        Args:
            original_material: The original material before transformation
            layer_thickness: Thickness of layers for fingerprint comparison
            threshold: Minimum improvement threshold to consider a swap detected

        Returns:
            Lattice: Corrected lattice with original orientation
        """
        swap_info = self.detect_swap_from_original(original_material, layer_thickness, threshold)

        if not swap_info.is_swapped:
            return self.material.lattice

        inverse_permutation = self._invert_permutation(swap_info.permutation)
        corrected_lattice = self._apply_permutation_to_lattice(self.material.lattice, inverse_permutation)
        return corrected_lattice

    def correct_material_to_match_original(
        self, original_material: MaterialWithBuildMetadata, layer_thickness: float = 1.0, threshold: float = 0.05
    ) -> MaterialWithBuildMetadata:
        """
        Correct the current material's lattice and basis to match the original material's orientation.

        This method detects any rotation/swap between the current and original material,
        and applies the inverse transformation to correct both the lattice vectors and
        the basis coordinates to preserve crystal coordinates.

        Args:
            original_material: The original material before transformation
            layer_thickness: Thickness of layers for fingerprint comparison
            threshold: Minimum improvement threshold to consider a swap detected

        Returns:
            MaterialWithBuildMetadata: Material with corrected lattice and transformed basis
        """
        corrected_material = self.material.clone()
        rotation_result = self.detect_rotation_from_original(original_material, layer_thickness, threshold)

        if not rotation_result.is_rotated:
            return corrected_material

        inverse_permutation = self._invert_permutation(
            self._rotation_matrix_to_permutation(rotation_result.rotation_matrix)
        )
        corrected_lattice = self._apply_permutation_to_lattice(self.material.lattice, inverse_permutation)

        corrected_material.set_lattice(corrected_lattice)
        corrected_material = rotate(corrected_material, axis=rotation_result.axis, angle=rotation_result.angle)
        corrected_material = wrap_to_unit_cell(corrected_material)

        return corrected_material
