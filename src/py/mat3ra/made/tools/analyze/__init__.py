from mat3ra.made.material import Material
from scipy.spatial.distance import pdist


class BaseMaterialAnalyzer:
    def __init__(self, material: Material):
        self.material = material

    @property
    def volume(self):
        return self.material.lattice.volume

    @property
    def atomic_density(self):
        return len(self.material.coordinates_array) / self.volume()

    def pairwise_distances(self):
        return pdist(self.material.coordinates_array)
