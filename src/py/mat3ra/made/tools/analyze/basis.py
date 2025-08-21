from typing import List, Tuple

from ..build_components.metadata import MaterialWithBuildMetadata
from . import BaseMaterialAnalyzer


class BasisMaterialAnalyzer(BaseMaterialAnalyzer):
    """Analyzer for material basis (atomic positions and elements)."""

    def get_layer_fingerprint(self, layer_thickness: float = 1.0) -> List[Tuple[Tuple[float, float], List[str]]]:
        """
        Create a fingerprint of the material by analyzing layers along the z-axis.

        Args:
            layer_thickness (float): Thickness of each layer in Angstroms

        Returns:
            List[Tuple[Tuple[float, float], List[str]]]: List of ((min_z, max_z), elements) for each layer
            Elements list contains sorted unique species in the layer, or [] if empty
        """
        material_cartesian = self.material.clone()
        material_cartesian.to_cartesian()

        coordinates = material_cartesian.basis.coordinates.values
        elements = material_cartesian.basis.elements.values

        z_coords = [coord[2] for coord in coordinates]
        min_z, max_z = min(z_coords), max(z_coords)

        layers = []
        current_z = min_z

        while current_z < max_z:
            layer_min = current_z
            layer_max = current_z + layer_thickness

            layer_elements = []
            for i, z in enumerate(z_coords):
                if layer_min <= z < layer_max:
                    layer_elements.append(elements[i])

            if layer_elements:
                unique_elements = sorted(list(set(layer_elements)))
            else:
                unique_elements = []  # Empty layer

            layers.append(((layer_min, layer_max), unique_elements))
            current_z += layer_thickness

        return layers

    def is_orientation_flipped(
        self, original_material: MaterialWithBuildMetadata, layer_thickness: float = 1.0
    ) -> bool:
        """
        Detect if the material orientation is flipped compared to the original.
        Uses Jaccard similarity to compare fingerprints in normal and flipped orientations.

        Args:
            original_material: The original material before primitivization
            layer_thickness: Thickness of layers for fingerprint comparison

        Returns:
            bool: True if orientation is flipped, False otherwise
        """
        original_analyzer = BasisMaterialAnalyzer(material=original_material)
        original_fingerprint = original_analyzer.get_layer_fingerprint(layer_thickness)
        current_fingerprint = self.get_layer_fingerprint(layer_thickness)

        # Compare fingerprints using similarity scores
        normal_score = self._compare_fingerprints(original_fingerprint, current_fingerprint)
        flipped_score = self._compare_fingerprints(original_fingerprint, list(reversed(current_fingerprint)))

        # If flipped orientation has significantly higher similarity, material is flipped
        # Use a threshold to determine if the difference is significant enough
        threshold = 0.1  # 10% difference threshold
        return flipped_score > normal_score + threshold

    def _compare_fingerprints(
        self, fp1: List[Tuple[Tuple[float, float], List[str]]], fp2: List[Tuple[Tuple[float, float], List[str]]]
    ) -> float:
        """
        Compare two fingerprints using Jaccard similarity.

        Args:
            fp1: First fingerprint
            fp2: Second fingerprint

        Returns:
            float: Average Jaccard similarity score (0.0 to 1.0)
        """
        if not fp1 or not fp2:
            return 0.0

        # Take the minimum length to avoid index errors
        min_length = min(len(fp1), len(fp2))
        if min_length == 0:
            return 0.0

        total_score = 0.0

        for i in range(min_length):
            _, elements1 = fp1[i]
            _, elements2 = fp2[i]

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
