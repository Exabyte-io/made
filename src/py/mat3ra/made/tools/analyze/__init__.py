from typing import Union

import numpy as np
from mat3ra.made.material import Material
from pydantic import BaseModel
from scipy.spatial.distance import pdist

from ..build_components.metadata.material_with_build_metadata import MaterialWithBuildMetadata
from .other import get_chemical_formula_empirical
from .utils import decorator_perform_operation_in_cartesian_coordinates


class BaseMaterialAnalyzer(BaseModel):
    material: Union[Material, MaterialWithBuildMetadata]

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

    @property
    def formula(self):
        return get_chemical_formula_empirical(self.material)
