import numpy as np
from mat3ra.made.material import Material
from scipy.spatial.distance import pdist


class BaseMaterialAnalyzer:
    def __init__(self, material: Material):
        self.material = material.clone()
        self.material.to_cartesian()

    @property
    def volume(self):
        return self.material.basis.cell.volume

    @property
    def atomic_density(self):
        return len(self.material.coordinates_array) / self.volume

    @property
    def pairwise_distances(self):
        return pdist(np.array(self.material.coordinates_array))
