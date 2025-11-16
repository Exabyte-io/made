from typing import ClassVar, Dict, List

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from pydantic import BaseModel, Field

from .layered_fingerprint_along_axis import LayeredFingerprintAlongAxis


class MaterialFingerprintAllAxes(BaseModel):
    """
    Complete fingerprint of a material across all three axes.
    """

    x_axis: LayeredFingerprintAlongAxis = Field(..., description="Fingerprint along x-axis")
    y_axis: LayeredFingerprintAlongAxis = Field(..., description="Fingerprint along y-axis")
    z_axis: LayeredFingerprintAlongAxis = Field(..., description="Fingerprint along z-axis")
    layer_thickness: float = Field(default=1.0, gt=0, description="Thickness of each layer in Angstroms")

    ALL_AXES: ClassVar[List[AxisEnum]] = [AxisEnum.x, AxisEnum.y, AxisEnum.z]

    def get_fingerprint_for_axis(self, axis: AxisEnum) -> LayeredFingerprintAlongAxis:
        axis_map = {
            AxisEnum.x: self.x_axis,
            AxisEnum.y: self.y_axis,
            AxisEnum.z: self.z_axis,
        }
        return axis_map[axis]

    def get_similarity_matrix(
        self, fingerprint_to_compare: "MaterialFingerprintAllAxes"
    ) -> Dict[str, Dict[str, float]]:
        """
        Calculate similarity matrix between this and another material fingerprint.

        Args:
            fingerprint_to_compare: Another MaterialFingerprintAllAxes to compare with

        Returns:
            Dict: Nested dictionary with similarity scores between all axis combinations
        """
        similarity_matrix = {}

        for self_axis in self.ALL_AXES:
            similarity_matrix[self_axis.value] = {}
            self_fp = self.get_fingerprint_for_axis(self_axis)

            for other_axis in self.ALL_AXES:
                other_fp = fingerprint_to_compare.get_fingerprint_for_axis(other_axis)
                similarity_matrix[self_axis.value][other_axis.value] = self_fp.get_similarity_score(other_fp)

        return similarity_matrix

    @staticmethod
    def reverse_axis_fingerprint(fingerprint: LayeredFingerprintAlongAxis) -> LayeredFingerprintAlongAxis:
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

