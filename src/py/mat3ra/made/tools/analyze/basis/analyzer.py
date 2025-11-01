from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.made.utils import AXIS_TO_INDEX_MAP

from ...build_components.metadata import MaterialWithBuildMetadata
from ...operations.core.unary import rotate
from .. import BaseMaterialAnalyzer
from .fingerprint import LayeredFingerprintAlongAxis, LayerFingerprint, MaterialFingerprint


class BasisMaterialAnalyzer(BaseMaterialAnalyzer):
    def get_layer_fingerprint(
        self, layer_thickness: float = 1.0, axis: AxisEnum = AxisEnum.z
    ) -> LayeredFingerprintAlongAxis:
        """
        Create a fingerprint of the material by analyzing layers along the specified axis.

        Args:
            layer_thickness: Thickness of each layer in Angstroms
            axis: Axis along which to compute the fingerprint

        Returns:
            LayeredFingerprintAlongAxis: Structured fingerprint with layer information
        """
        material_cartesian = self.material.clone()
        material_cartesian.to_cartesian()

        coordinates = material_cartesian.basis.coordinates.values
        elements = material_cartesian.basis.elements.values

        axis_index = AXIS_TO_INDEX_MAP[axis.value]
        axis_coords = [coord[axis_index] for coord in coordinates]
        min_coord, max_coord = min(axis_coords), max(axis_coords)

        fingerprint = LayeredFingerprintAlongAxis(axis=axis, layer_thickness=layer_thickness)

        current_coord = min_coord
        while current_coord < max_coord:
            layer_min = current_coord
            layer_max = current_coord + layer_thickness

            layer_elements = []
            for i, coord_val in enumerate(axis_coords):
                if layer_min <= coord_val < layer_max:
                    layer_elements.append(elements[i])

            unique_elements = sorted(list(set(layer_elements))) if layer_elements else []
            layer = LayerFingerprint(min_coord=layer_min, max_coord=layer_max, elements=unique_elements)
            fingerprint.layers.append(layer)

            current_coord += layer_thickness

        return fingerprint

    def get_material_fingerprint(self, layer_thickness: float = 1.0) -> MaterialFingerprint:
        """
        Create a complete fingerprint of the material across all three axes.

        Args:
            layer_thickness: Thickness of each layer in Angstroms

        Returns:
            MaterialFingerprint: Complete fingerprint with layer information for all axes
        """
        x_fingerprint = self.get_layer_fingerprint(layer_thickness, AxisEnum.x)
        y_fingerprint = self.get_layer_fingerprint(layer_thickness, AxisEnum.y)
        z_fingerprint = self.get_layer_fingerprint(layer_thickness, AxisEnum.z)

        return MaterialFingerprint(
            x_axis=x_fingerprint,
            y_axis=y_fingerprint,
            z_axis=z_fingerprint,
            layer_thickness=layer_thickness
        )

    def detect_rotation_from_original(
        self, original_material: MaterialWithBuildMetadata, layer_thickness: float = 1.0, threshold: float = 0.1
    ) -> dict:
        """
        Detect rotation of the current material compared to the original material.
        
        Args:
            original_material: The original material before transformation
            layer_thickness: Thickness of layers for fingerprint comparison
            threshold: Minimum improvement threshold to consider a rotation detected
            
        Returns:
            dict: Rotation detection results with rotation type, axis, and confidence
        """
        original_analyzer = BasisMaterialAnalyzer(material=original_material)
        original_fingerprint = original_analyzer.get_material_fingerprint(layer_thickness)
        current_fingerprint = self.get_material_fingerprint(layer_thickness)
        
        return original_fingerprint.detect_rotation(current_fingerprint, threshold)

    def apply_corrective_rotation(
        self, original_material: MaterialWithBuildMetadata, layer_thickness: float = 1.0
    ) -> MaterialWithBuildMetadata:
        """
        Apply corrective rotation to align the material with the original orientation.
        
        Args:
            original_material: The original material before transformation
            layer_thickness: Thickness of layers for fingerprint comparison
            
        Returns:
            MaterialWithBuildMetadata: Material with corrective rotation applied
        """
        rotation_info = self.detect_rotation_from_original(original_material, layer_thickness)
        
        if not rotation_info['is_rotated']:
            return self.material
            
        # Apply the inverse rotation to correct the orientation
        rotation_axis = rotation_info['rotation_axis']
        rotation_angle = -rotation_info['rotation_angle']  # Inverse rotation
        
        corrected_material = rotate(
            self.material, 
            axis=rotation_axis.tolist(), 
            angle=rotation_angle, 
            rotate_cell=False
        )
                
        return corrected_material

