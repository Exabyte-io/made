from typing import List

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
