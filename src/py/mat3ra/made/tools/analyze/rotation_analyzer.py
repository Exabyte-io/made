from typing import List, Optional, Union

import numpy as np
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.made.material import Material
from pydantic import BaseModel, Field
from scipy.spatial.transform import Rotation

from .basis.analyzer import BasisMaterialAnalyzer
from .fingerprint import MaterialFingerprintAllAxes
from ..build_components.metadata.material_with_build_metadata import MaterialWithBuildMetadata


class RotationParameters(BaseModel):
    model_config = {"arbitrary_types_allowed": True}

    rotation_matrix: np.ndarray = Field(..., description="3x3 rotation matrix")
    angle: float = Field(..., description="Rotation angle in degrees")
    axis: List[int] = Field(..., description="Rotation axis as [x, y, z] for rotate() function")


class RotationDetectionResult(RotationParameters):
    is_rotated: bool = Field(..., description="Whether a rotation was detected")
    confidence: float = Field(..., description="Confidence score of the detection (0.0 to 1.0)")


class MaterialRotationAnalyzer(BaseModel):
    material: Union[Material, MaterialWithBuildMetadata]
    rotation_axis_significance_threshold: float = 0.5

    def detect_rotation_from_original(
        self, original_material: MaterialWithBuildMetadata, layer_thickness: float = 1.0, threshold: float = 0.1
    ) -> RotationDetectionResult:
        """
        Detect rotation of the current material compared to the original material.

        Args:
            original_material: The original material before transformation
            layer_thickness: Thickness of layers for fingerprint comparison
            threshold: Minimum improvement threshold to consider a rotation detected

        Returns:
            RotationDetectionResult: Rotation detection results with rotation type, axis, and confidence
        """
        original_analyzer = BasisMaterialAnalyzer(material=original_material)
        original_fingerprint = original_analyzer.get_material_fingerprint(layer_thickness)

        current_analyzer = BasisMaterialAnalyzer(material=self.material)
        current_fingerprint = current_analyzer.get_material_fingerprint(layer_thickness)

        return self._detect_rotation_between_fingerprints(original_fingerprint, current_fingerprint, threshold)

    def _detect_rotation_between_fingerprints(
        self,
        original_fingerprint: MaterialFingerprintAllAxes,
        current_fingerprint: MaterialFingerprintAllAxes,
        threshold: float = 0.1,
    ) -> RotationDetectionResult:
        """
        Detect rotation between two material fingerprints.

        Args:
            original_fingerprint: Original material fingerprint
            current_fingerprint: Current material fingerprint
            threshold: Minimum improvement threshold to consider a rotation detected

        Returns:
            RotationDetectionResult: Rotation detection results
        """
        direct_score = self._calculate_alignment_score(
            original_fingerprint, current_fingerprint, rotation_matrix=np.eye(3)
        )

        rotation_candidates = self._generate_rotation_candidates()

        best_score = direct_score
        best_result = RotationDetectionResult(
            is_rotated=False,
            rotation_matrix=np.eye(3),
            angle=0.0,
            axis=[0, 0, 0],
            confidence=direct_score,
        )

        for rotation_params in rotation_candidates:
            score = self._calculate_alignment_score(
                original_fingerprint, current_fingerprint, rotation_params.rotation_matrix
            )

            if score > best_score + threshold:
                best_score = score
                best_result = RotationDetectionResult(
                    is_rotated=True,
                    rotation_matrix=rotation_params.rotation_matrix,
                    angle=rotation_params.angle,
                    axis=rotation_params.axis,
                    confidence=score,
                )

        return best_result

    def _calculate_alignment_score(
        self,
        original_fingerprint: MaterialFingerprintAllAxes,
        current_fingerprint: MaterialFingerprintAllAxes,
        rotation_matrix: np.ndarray,
    ) -> float:
        """
        Calculate alignment score between fingerprints with rotation.

        Args:
            original_fingerprint: Original material fingerprint
            current_fingerprint: Current material fingerprint
            rotation_matrix: Rotation matrix to apply before comparison (use np.eye(3) for no rotation)

        Returns:
            float: Alignment score (0.0 to 1.0)
        """
        total_score = 0.0
        count = 0

        for axis_index, axis_name in enumerate(MaterialFingerprintAllAxes.ALL_AXES):
            axis_score = self._calculate_axis_alignment_score(
                axis_index,
                axis_name,
                original_fingerprint,
                current_fingerprint,
                rotation_matrix,
            )

            if axis_score is not None:
                total_score += axis_score
                count += 1

        return total_score / count if count > 0 else 0.0

    def _calculate_axis_alignment_score(
        self,
        axis_index: int,
        axis_name: AxisEnum,
        original_fingerprint: MaterialFingerprintAllAxes,
        current_fingerprint: MaterialFingerprintAllAxes,
        rotation_matrix: np.ndarray,
    ) -> Optional[float]:
        row = rotation_matrix[axis_index, :]
        max_idx = int(np.argmax(np.abs(row)))

        if np.abs(row[max_idx]) <= self.rotation_axis_significance_threshold:
            return None

        reference_fingerprint = original_fingerprint.get_fingerprint_for_axis(axis_name)
        target_axis = MaterialFingerprintAllAxes.ALL_AXES[max_idx]
        target_fingerprint = current_fingerprint.get_fingerprint_for_axis(target_axis)

        if row[max_idx] < 0:
            target_fingerprint = MaterialFingerprintAllAxes.reverse_axis_fingerprint(target_fingerprint)

        return reference_fingerprint.get_similarity_score(target_fingerprint)

    def _generate_rotation_candidates(self) -> List[RotationParameters]:
        """
        Generate common rotation matrices for testing.

        Returns:
            List of RotationParameters with rotation matrix, angle, and axis
        """
        candidates = []
        angles = [90, -90, 180]
        axes = np.eye(3, dtype=int)

        for axis_vector in axes:
            for angle in angles:
                rotation_matrix = Rotation.from_rotvec(axis_vector * np.radians(angle)).as_matrix()
                candidates.append(
                    RotationParameters(rotation_matrix=rotation_matrix, angle=angle, axis=axis_vector.tolist())
                )
        return candidates
