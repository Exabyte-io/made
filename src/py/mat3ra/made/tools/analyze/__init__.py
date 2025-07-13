import numpy as np
from pydantic import BaseModel
from scipy.spatial.distance import pdist

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.utils import decorator_perform_operation_in_cartesian_coordinates
from mat3ra.made.tools.build import MaterialWithBuildMetadata


class BaseMaterialAnalyzer(BaseModel):
    material: Material

    @property
    def volume(self):
        return self.material.basis.cell.volume

    @property
    def atomic_density(self):
        return len(self.material.coordinates_array) / self.volume

    @property
    @decorator_perform_operation_in_cartesian_coordinates
    def pairwise_distances(self):
        return pdist(np.array(self.material.coordinates_array))
