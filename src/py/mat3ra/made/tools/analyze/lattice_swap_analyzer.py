import math
from typing import Union

import numpy as np
from mat3ra.code.array_with_ids import ArrayWithIds
from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from pydantic import BaseModel, Field

from mat3ra.made.material import Material
from .basis import BasisMaterialAnalyzer
from .enums import POSSIBLE_TRANSFORMATION_MATRICES
from ..build_components.metadata.material_with_build_metadata import MaterialWithBuildMetadata
from ..operations.reusable.unary import transform_material_by_matrix


class LatticeSwapDetectionResult(BaseModel):
    model_config = {"arbitrary_types_allowed": True}
    is_swapped: bool = Field(..., description="Whether a lattice swap was detected")
    permutation: Matrix3x3Schema = Field(..., description="Transformation matrix representing the swap")
    confidence: float = Field(..., description="Confidence score (0.0 to 1.0)")


class MaterialLatticeSwapAnalyzer(BaseModel):
    """
    Analyzer to detect lattice vector swaps/permutations between materials.

    This detects when lattice vectors have been reoriented (e.g., a->a, b->c, c->-b)
    rather than the basis being rotated within the lattice.
    """

    material: Union[Material, MaterialWithBuildMetadata]
    tolerance: float = 0.01

    def _compute_transformation_score(self, matrix: np.ndarray, target_fingerprint) -> float:
        transformed_material = transform_material_by_matrix(self.material, matrix)
        new_analyzer = BasisMaterialAnalyzer(material=transformed_material)
        new_fingerprint = new_analyzer.get_layer_fingerprint(target_fingerprint.layer_thickness)
        return target_fingerprint.get_similarity_score(new_fingerprint)

    def _create_detection_result(self, matrix: np.ndarray, score: float) -> LatticeSwapDetectionResult:
        is_identity = np.allclose(matrix, np.eye(3), atol=self.tolerance)
        return LatticeSwapDetectionResult(
            is_swapped=not is_identity,
            permutation=Matrix3x3Schema(root=matrix.tolist()),
            confidence=score,
        )

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
        target_analyzer = BasisMaterialAnalyzer(material=original_material)
        target_fingerprint = target_analyzer.get_layer_fingerprint(layer_thickness)
        possible_matrices = list(map(np.array, POSSIBLE_TRANSFORMATION_MATRICES))

        best_match = None
        max_score = 0

        for matrix in possible_matrices:
            score = self._compute_transformation_score(matrix, target_fingerprint)
            if score > max_score:
                best_match = self._create_detection_result(matrix, score)
                max_score = score

        if best_match and best_match.is_swapped and best_match.confidence >= threshold:
            return best_match
        # If no swap detected, return identity result
        return self._create_detection_result(np.eye(3), 1.0)

    def get_corrected_material(
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
            matrix_list = [list(row.root) if hasattr(row, "root") else row for row in swap_info.permutation.root]
            matrix_array = np.array(matrix_list)
            return self.canonicalize_in_plane_material(transform_material_by_matrix(self.material, matrix_array))
        return self.canonicalize_in_plane_material(self.material)

    @classmethod
    def canonicalize_in_plane_material(
        cls, material: MaterialWithBuildMetadata, tolerance: float = 1e-8
    ) -> MaterialWithBuildMetadata:
        if not cls._should_canonicalize_in_plane(material, tolerance):
            return material

        canonical_material = material.clone()
        original_is_cartesian = canonical_material.basis.is_in_cartesian_units
        canonical_material.to_crystal()

        coordinates = np.array(canonical_material.coordinates_array, dtype=float)
        if len(coordinates) == 0:
            return canonical_material

        elements = list(canonical_material.basis.elements.values)
        labels = list(canonical_material.basis.labels.values)
        constraints = list(canonical_material.basis.constraints.values)

        best_representation = None
        best_coordinates = None
        best_order = None

        for operation in cls._get_in_plane_symmetry_operations(material, tolerance):
            transformed_coordinates = coordinates.copy()
            for index, coordinate in enumerate(coordinates):
                x_value, y_value = operation(coordinate[0], coordinate[1])
                transformed_coordinates[index, 0] = x_value
                transformed_coordinates[index, 1] = y_value

            x_candidates = sorted(
                {
                    cls._normalize_fractional_coordinate(coord[0], tolerance)
                    for coord in transformed_coordinates
                }
            )
            y_candidates = sorted(
                {
                    cls._normalize_fractional_coordinate(coord[1], tolerance)
                    for coord in transformed_coordinates
                }
            )

            for x_shift in x_candidates:
                for y_shift in y_candidates:
                    shifted_coordinates = transformed_coordinates.copy()
                    shifted_coordinates[:, 0] = [
                        cls._normalize_fractional_coordinate(value - x_shift, tolerance)
                        for value in shifted_coordinates[:, 0]
                    ]
                    shifted_coordinates[:, 1] = [
                        cls._normalize_fractional_coordinate(value - y_shift, tolerance)
                        for value in shifted_coordinates[:, 1]
                    ]

                    order = sorted(
                        range(len(shifted_coordinates)),
                        key=lambda index: cls._canonical_atom_key(
                            elements[index],
                            labels[index] if labels else None,
                            shifted_coordinates[index],
                            tolerance,
                        ),
                    )
                    representation = tuple(
                        cls._canonical_atom_key(
                            elements[index],
                            labels[index] if labels else None,
                            shifted_coordinates[index],
                            tolerance,
                        )
                        for index in order
                    )

                    if best_representation is None or representation < best_representation:
                        best_representation = representation
                        best_coordinates = shifted_coordinates
                        best_order = order

        if best_coordinates is None or best_order is None:
            return material

        canonical_material.basis.coordinates.values = best_coordinates[best_order].tolist()
        canonical_material.basis.elements = ArrayWithIds.from_values(
            [elements[index] for index in best_order]
        )
        if labels:
            canonical_material.basis.labels = ArrayWithIds.from_values(
                [labels[index] for index in best_order]
            )
        if constraints and len(constraints) == len(best_order):
            canonical_material.basis.constraints = ArrayWithIds.from_values(
                [constraints[index] for index in best_order]
            )

        if original_is_cartesian:
            canonical_material.to_cartesian()

        return canonical_material

    @staticmethod
    def _should_canonicalize_in_plane(material: MaterialWithBuildMetadata, tolerance: float) -> bool:
        return (
            math.isclose(material.lattice.alpha, 90.0, abs_tol=tolerance)
            and math.isclose(material.lattice.beta, 90.0, abs_tol=tolerance)
            and math.isclose(material.lattice.gamma, 90.0, abs_tol=tolerance)
            and math.isclose(material.lattice.a, material.lattice.b, rel_tol=0.0, abs_tol=tolerance)
            and material.lattice.c > max(material.lattice.a, material.lattice.b)
        )

    @staticmethod
    def _get_in_plane_symmetry_operations(material: MaterialWithBuildMetadata, tolerance: float):
        operations = [
            lambda x_value, y_value: (x_value, y_value),
            lambda x_value, y_value: (-x_value, y_value),
            lambda x_value, y_value: (x_value, -y_value),
            lambda x_value, y_value: (-x_value, -y_value),
        ]

        if math.isclose(material.lattice.a, material.lattice.b, rel_tol=0.0, abs_tol=tolerance):
            operations.extend(
                [
                    lambda x_value, y_value: (y_value, x_value),
                    lambda x_value, y_value: (-y_value, x_value),
                    lambda x_value, y_value: (y_value, -x_value),
                    lambda x_value, y_value: (-y_value, -x_value),
                ]
            )

        return operations

    @staticmethod
    def _normalize_fractional_coordinate(value: float, tolerance: float) -> float:
        normalized_value = float(value % 1.0)
        if abs(normalized_value) < tolerance or abs(normalized_value - 1.0) < tolerance:
            return 0.0
        return normalized_value

    @classmethod
    def _canonical_atom_key(cls, element, label, coordinate, tolerance: float):
        rounded_coordinate = tuple(
            round(cls._normalize_fractional_coordinate(value, tolerance), 8)
            for value in coordinate
        )
        label_key = "" if label is None else str(label)
        return (element, label_key, rounded_coordinate[2], rounded_coordinate[0], rounded_coordinate[1])
