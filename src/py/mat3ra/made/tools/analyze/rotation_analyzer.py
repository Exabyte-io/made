from typing import List, Optional, Tuple, Union

import numpy as np
from mat3ra.made.material import Material
from mat3ra.made.utils import AXIS_TO_INDEX_MAP
from pydantic import BaseModel, Field
from scipy.spatial.transform import Rotation

from ..build_components.metadata.material_with_build_metadata import MaterialWithBuildMetadata
from .basis.analyzer import BasisMaterialAnalyzer
from .fingerprint import LayeredFingerprintAlongAxis, MaterialFingerprintAllAxes


class RotationDetectionResult(BaseModel):
    model_config = {"arbitrary_types_allowed": True}

    is_rotated: bool = Field(..., description="Whether a rotation was detected")
    rotation_matrix: Optional[np.ndarray] = Field(
        None, description="3x3 rotation matrix if rotation detected, None otherwise"
    )
    rotation_angle: Optional[float] = Field(None, description="Rotation angle in degrees if detected, None otherwise")
    rotation_axis: Optional[List[int]] = Field(
        None, description="Rotation axis as [x, y, z] for rotate() function if detected, None otherwise"
    )
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
            rotation_matrix=None,
            rotation_angle=None,
            rotation_axis=None,
            confidence=direct_score,
        )

        for rotation_matrix, angle, axis in rotation_candidates:
            score = self._calculate_alignment_score(original_fingerprint, current_fingerprint, rotation_matrix)

            if score > best_score + threshold:
                best_score = score
                best_result = RotationDetectionResult(
                    is_rotated=True,
                    rotation_matrix=rotation_matrix,
                    rotation_angle=angle,
                    rotation_axis=axis,
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

        for i, self_axis in enumerate(MaterialFingerprintAllAxes.ALL_AXES):
            self_fp = original_fingerprint.get_fingerprint_for_axis(self_axis)

            row = rotation_matrix[i, :]
            max_idx = np.argmax(np.abs(row))

            if np.abs(row[max_idx]) > self.rotation_axis_significance_threshold:
                other_axis = MaterialFingerprintAllAxes.ALL_AXES[max_idx]
                other_fp = current_fingerprint.get_fingerprint_for_axis(other_axis)

                if row[max_idx] < 0:
                    other_fp = self._reverse_axis_fingerprint(other_fp)

                total_score += self_fp.get_similarity_score(other_fp)
                count += 1

        return total_score / count if count > 0 else 0.0

    def _reverse_axis_fingerprint(self, fingerprint: LayeredFingerprintAlongAxis) -> LayeredFingerprintAlongAxis:
        """
        Reverse a fingerprint along its axis (for 180-degree rotations).

        Args:
            fingerprint: Fingerprint to reverse

        Returns:
            LayeredFingerprintAlongAxis: Reversed fingerprint
        """
        reversed_layers = list(reversed(fingerprint.layers))
        return LayeredFingerprintAlongAxis(
            layers=reversed_layers, axis=fingerprint.axis, layer_thickness=fingerprint.layer_thickness
        )

    def _generate_rotation_candidates(self) -> List[Tuple[np.ndarray, float, List[int]]]:
        """
        Generate common rotation matrices for testing.

        Returns:
            List of tuples (rotation_matrix, angle_degrees, axis_list)
        """
        candidates = []

        for axis_enum in MaterialFingerprintAllAxes.ALL_AXES:
            axis_idx = AXIS_TO_INDEX_MAP[axis_enum.value]
            for angle in [90, -90, 180]:
                axis_vector = np.zeros(3)
                axis_vector[axis_idx] = 1.0

                axis_list = [int(axis_vector[0]), int(axis_vector[1]), int(axis_vector[2])]

                rotation_matrix = Rotation.from_rotvec(axis_vector * np.radians(angle)).as_matrix()
                candidates.append((rotation_matrix, angle, axis_list))

        return candidates
