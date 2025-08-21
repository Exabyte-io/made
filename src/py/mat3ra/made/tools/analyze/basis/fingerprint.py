from typing import List

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from pydantic import BaseModel, Field


class LayerFingerprint(BaseModel):
    """
    Represents a single layer in the material fingerprint.
    """

    min_coord: float = Field(..., description="Minimum coordinate value for the layer")
    max_coord: float = Field(..., description="Maximum coordinate value for the layer")
    elements: List[str] = Field(default_factory=list, description="Sorted unique chemical elements in the layer")


class LayeredFingerprintAlongAxis(BaseModel):
    """
    Material fingerprint along a specific axis.
    """

    layers: List[LayerFingerprint] = Field(default_factory=list, description="List of layer fingerprints")
    axis: AxisEnum = Field(default=AxisEnum.z, description="Axis along which the fingerprint is computed")
    layer_thickness: float = Field(default=1.0, gt=0, description="Thickness of each layer in Angstroms")

    class Config:
        """Pydantic configuration."""

        use_enum_values = True
        validate_assignment = True

    def get_non_empty_layers(self) -> List[LayerFingerprint]:
        """Get only layers that contain chemical elements."""
        return [layer for layer in self.layers if layer.elements]

    def get_element_sequence(self) -> List[List[str]]:
        """Get element lists for all layers (including empty)."""
        return [layer.elements for layer in self.layers]

    def get_non_empty_element_sequence(self) -> List[List[str]]:
        """Get element lists for non-empty layers only."""
        return [layer.elements for layer in self.layers if layer.elements]
