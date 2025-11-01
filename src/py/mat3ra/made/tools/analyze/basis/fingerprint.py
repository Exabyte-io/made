from typing import ClassVar, Dict, List, Tuple

import numpy as np
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from pydantic import BaseModel, Field


class LayerFingerprint(BaseModel):
    min_coord: float = Field(..., description="Minimum coordinate value for the layer")
    max_coord: float = Field(..., description="Maximum coordinate value for the layer")
    elements: List[str] = Field(default_factory=list, description="Sorted unique chemical elements in the layer")


class LayeredFingerprintAlongAxis(BaseModel):
    layers: List[LayerFingerprint] = Field(default_factory=list, description="List of layer fingerprints")
    axis: AxisEnum = Field(default=AxisEnum.z, description="Axis along which the fingerprint is computed")
    layer_thickness: float = Field(default=1.0, gt=0, description="Thickness of each layer in Angstroms")

    def get_non_empty_layers(self) -> List[LayerFingerprint]:
        return [layer for layer in self.layers if layer.elements]

    def get_element_sequence(self) -> List[List[str]]:
        return [layer.elements for layer in self.layers]

    def get_non_empty_element_sequence(self) -> List[List[str]]:
        return [layer.elements for layer in self.layers if layer.elements]

    def get_similarity_score(self, other: "LayeredFingerprintAlongAxis") -> float:
        """
        Calculate Jaccard similarity score between this and another fingerprint.

        Args:
            other: Another LayeredFingerprintAlongAxis to compare with

        Returns:
            float: Average Jaccard similarity score (0.0 to 1.0)
        """
        if not self.layers or not other.layers:
            return 0.0

        min_length = min(len(self.layers), len(other.layers))
        if min_length == 0:
            return 0.0

        total_score = 0.0

        for i in range(min_length):
            elements1 = self.layers[i].elements
            elements2 = other.layers[i].elements

            # Handle empty layers
            if len(elements1) == 0 and len(elements2) == 0:
                layer_score = 1.0
            elif len(elements1) == 0 or len(elements2) == 0:
                layer_score = 0.0
            else:
                # Calculate Jaccard similarity
                set1 = set(elements1)
                set2 = set(elements2)
                intersection = len(set1.intersection(set2))
                union = len(set1.union(set2))
                layer_score = intersection / union if union > 0 else 0.0

            total_score += layer_score

        return total_score / min_length


class MaterialFingerprint(BaseModel):
    """
    Complete fingerprint of a material across all three axes.
    """
    x_axis: LayeredFingerprintAlongAxis = Field(..., description="Fingerprint along x-axis")
    y_axis: LayeredFingerprintAlongAxis = Field(..., description="Fingerprint along y-axis")
    z_axis: LayeredFingerprintAlongAxis = Field(..., description="Fingerprint along z-axis")
    layer_thickness: float = Field(default=1.0, gt=0, description="Thickness of each layer in Angstroms")
    
    ROTATION_AXIS_SIGNIFICANCE_THRESHOLD: ClassVar[float] = 0.5  # Minimum component magnitude to consider axis mapping significant

    def get_fingerprint_for_axis(self, axis: AxisEnum) -> LayeredFingerprintAlongAxis:
        axis_map = {
            AxisEnum.x: self.x_axis,
            AxisEnum.y: self.y_axis,
            AxisEnum.z: self.z_axis,
        }
        return axis_map[axis]

    def get_all_axis_fingerprints(self) -> Dict[AxisEnum, LayeredFingerprintAlongAxis]:
        return {
            AxisEnum.x: self.x_axis,
            AxisEnum.y: self.y_axis,
            AxisEnum.z: self.z_axis,
        }

    def get_similarity_matrix(self, other: "MaterialFingerprint") -> Dict[str, Dict[str, float]]:
        """
        Calculate similarity matrix between this and another material fingerprint.
        
        Args:
            other: Another MaterialFingerprint to compare with
            
        Returns:
            Dict: Nested dictionary with similarity scores between all axis combinations
        """
        axes = [AxisEnum.x, AxisEnum.y, AxisEnum.z]
        similarity_matrix = {}
        
        for self_axis in axes:
            similarity_matrix[self_axis.value] = {}
            self_fp = self.get_fingerprint_for_axis(self_axis)
            
            for other_axis in axes:
                other_fp = other.get_fingerprint_for_axis(other_axis)
                similarity_matrix[self_axis.value][other_axis.value] = self_fp.get_similarity_score(other_fp)
                
        return similarity_matrix

    def detect_rotation(self, other: "MaterialFingerprint", threshold: float = 0.1) -> Dict[str, any]:
        """
        Detect rotation between this and another material fingerprint.
        
        Args:
            other: Another MaterialFingerprint to compare with
            threshold: Minimum improvement threshold to consider a rotation detected
            
        Returns:
            Dict containing rotation information: {
                'is_rotated': bool,
                'rotation_matrix': np.ndarray or None,  # 3x3 rotation matrix
                'rotation_angle': float or None,        # rotation angle in degrees
                'rotation_axis': np.ndarray or None,    # rotation axis vector
                'confidence': float                     # confidence score of the detection
            }
        """
        direct_score = self._calculate_alignment_score(other, identity=True)
        
        rotation_candidates = self._generate_rotation_candidates()
        
        best_score = direct_score
        best_info = {
            'is_rotated': False,
            'rotation_matrix': None,
            'rotation_angle': None,
            'rotation_axis': None,
            'confidence': direct_score
        }
        
        for rotation_matrix, angle, axis in rotation_candidates:
            score = self._calculate_alignment_score(other, rotation_matrix=rotation_matrix)
            
            if score > best_score + threshold:
                best_score = score
                best_info = {
                    'is_rotated': True,
                    'rotation_matrix': rotation_matrix,
                    'rotation_angle': angle,
                    'rotation_axis': axis,
                    'confidence': score
                }
        
        return best_info

    def _calculate_alignment_score(self, other: "MaterialFingerprint", identity: bool = False, rotation_matrix: np.ndarray = None) -> float:
        """
        Calculate alignment score between fingerprints, optionally with rotation.
        
        Args:
            other: Another MaterialFingerprint to compare with
            identity: If True, compare directly without rotation
            rotation_matrix: Rotation matrix to apply before comparison
            
        Returns:
            float: Alignment score (0.0 to 1.0)
        """
        axes = [AxisEnum.x, AxisEnum.y, AxisEnum.z]
        
        if identity:
            total_score = 0.0
            for axis in axes:
                self_fp = self.get_fingerprint_for_axis(axis)
                other_fp = other.get_fingerprint_for_axis(axis)
                total_score += self_fp.get_similarity_score(other_fp)
            return total_score / len(axes)
        
        elif rotation_matrix is not None:
            total_score = 0.0
            count = 0
            
            for i, self_axis in enumerate(axes):
                self_fp = self.get_fingerprint_for_axis(self_axis)
                
                row = rotation_matrix[i, :]
                max_idx = np.argmax(np.abs(row))
                
                if np.abs(row[max_idx]) > self.ROTATION_AXIS_SIGNIFICANCE_THRESHOLD:
                    other_axis = axes[max_idx]
                    other_fp = other.get_fingerprint_for_axis(other_axis)
                    
                    # If the component is negative, we need to reverse the fingerprint
                    if row[max_idx] < 0:
                        other_fp = self._reverse_axis_fingerprint(other_fp)
                    
                    total_score += self_fp.get_similarity_score(other_fp)
                    count += 1
            
            return total_score / count if count > 0 else 0.0
        
        return 0.0

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
            layers=reversed_layers, 
            axis=fingerprint.axis, 
            layer_thickness=fingerprint.layer_thickness
        )

    def _generate_rotation_candidates(self) -> List[Tuple[np.ndarray, float, np.ndarray]]:
        """
        Generate common rotation matrices for testing.
        
        Returns:
            List of tuples (rotation_matrix, angle_degrees, axis_vector)
        """
        candidates = []
        
        for axis_idx, axis_name in enumerate(['x', 'y', 'z']):
            for angle in [90, -90, 180]:
                axis_vector = np.zeros(3)
                axis_vector[axis_idx] = 1.0
                
                rotation_matrix = self._rotation_matrix_from_axis_angle(axis_vector, np.radians(angle))
                candidates.append((rotation_matrix, angle, axis_vector))
        
        return candidates

    def _rotation_matrix_from_axis_angle(self, axis: np.ndarray, angle: float) -> np.ndarray:
        """
        Create rotation matrix from axis and angle using Rodrigues' formula.
        
        Args:
            axis: Unit vector representing rotation axis
            angle: Rotation angle in radians
            
        Returns:
            3x3 rotation matrix
        """
        axis = axis / np.linalg.norm(axis)
        cos_angle = np.cos(angle)
        sin_angle = np.sin(angle)
        
        K = np.array([
            [0, -axis[2], axis[1]],
            [axis[2], 0, -axis[0]],
            [-axis[1], axis[0], 0]
        ])
        
        R = np.eye(3) + sin_angle * K + (1 - cos_angle) * (K @ K)
        return R
