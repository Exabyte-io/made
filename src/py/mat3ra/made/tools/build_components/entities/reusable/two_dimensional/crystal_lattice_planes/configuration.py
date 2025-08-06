from typing import Union

from mat3ra.esse.models.materials_category_components.entities.reusable.two_dimensional.crystal_lattice_planes import (
    CrystalLatticePlanesSchema,
)
from mat3ra.made.material import Material
from mat3ra.made.tools.build_components.entities.reusable.base_builder import BaseConfigurationPydantic

from ......build_components import MaterialWithBuildMetadata
from ....auxiliary.two_dimensional.miller_indices import MillerIndices


class CrystalLatticePlanesConfiguration(CrystalLatticePlanesSchema, BaseConfigurationPydantic):
    crystal: Union[Material, MaterialWithBuildMetadata]

    @property
    def in_plane_vectors(self):
        # Two vectors in the plane of the Miller indices
        return self.crystal.lattice.vector_arrays[:2, :2]

    @property
    def miller_indices_as_string(self) -> str:
        miller_indices_cls_instance = MillerIndices(root=self.miller_indices)
        return str(miller_indices_cls_instance)
