import numpy as np
from mat3ra.made.material import Material
from pydantic import BaseModel

from . import BaseMaterialAnalyzer


class RadialDistributionFunction(BaseModel):
    rdf: np.ndarray
    bin_centers: np.ndarray

    class Config:
        arbitrary_types_allowed = True

    @classmethod
    def from_material(cls, material: Material, cutoff: float = 10.0, bin_size: float = 0.1):
        analyzer = BaseMaterialAnalyzer(material)
        distances = analyzer.pairwise_distances
        density = analyzer.atomic_density
        distances = distances[distances <= cutoff]

        # Bin distances into a histogram
        bins = np.arange(0, cutoff + bin_size, bin_size)  # Bin edges
        hist, bin_edges = np.histogram(distances, bins=bins, density=False)

        # Convert to radial distribution function
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        shell_volumes = (
            (4 / 3) * np.pi * (np.power(bin_edges[1:], 3) - np.power(bin_edges[:-1], 3))
        )  # Volume of spherical shells

        rdf = hist / (shell_volumes * density)

        return cls(rdf=rdf, bin_centers=bin_centers)

    @property
    def first_peak_index(self):
        return np.argmax(self.rdf[1:]) + 1

    @property
    def first_peak_value(self):
        return self.rdf[self.first_peak_index]

    @property
    def first_peak_width(self):
        half_max = 0.5 * self.first_peak_value

        # Find left boundary
        left_index = self.first_peak_index
        while left_index > 0 and self.rdf[left_index] > half_max:
            left_index -= 1

        # Find right boundary
        right_index = self.first_peak_index
        while right_index < len(self.rdf) - 1 and self.rdf[right_index] > half_max:
            right_index += 1

        # Compute the width
        left_boundary = self.bin_centers[left_index]
        right_boundary = self.bin_centers[right_index]
        first_peak_width = right_boundary - left_boundary

        return float(first_peak_width)

    @property
    def first_peak_distance(self):
        return self.bin_centers[self.first_peak_index]

    def is_within_first_peak(self, distance: float, tolerance: float = 0.1) -> bool:
        return (
            self.first_peak_distance - 0.5 * self.first_peak_width - tolerance
            < distance
            < self.first_peak_distance + 0.5 * self.first_peak_width + tolerance
        )
