from typing import Any

from mat3ra.esse.models.apse.materials.builders.slab.pymatgen.parameters import (
    PymatgenSlabGeneratorParametersSchema,
)
from pydantic import BaseModel

from mat3ra.made.material import Material
from .configuration import (
    SlabConfiguration,
    AtomicLayersUniqueRepeatedConfiguration,
    VacuumConfiguration,
    AtomicLayersUniqueRepeatedBuilder,
)
from .termination import Termination
from ..stack.builders import StackBuilder2Components
from ..vacuum.builders import VacuumBuilder
from ...build import BaseBuilder
from ...operations.core.unary import stack, supercell


class SlabSelectorParameters(BaseModel):
    termination: Termination


class PymatgenSlabGeneratorParameters(PymatgenSlabGeneratorParametersSchema):
    # Parameters described in https://github.com/materialsproject/pymatgen/blob/585bb673c4aa222669c4b0d72ffeec3dbf092630/pymatgen/core/surface.py#L1187
    pass


class SlabBuilderParameters(PymatgenSlabGeneratorParameters):
    make_primitive: bool = False
    use_orthogonal_c: bool = True


class SlabBuilder(StackBuilder2Components):

    def generate(self, configuration: SlabConfiguration) -> Material:
        stacked_materials = super().generate(configuration)
        supercell_slab = supercell(stacked_materials, configuration.xy_supercell_matrix)
        return supercell_slab

    def get_material(self, configuration: SlabConfiguration) -> Material:
        return self.generate(configuration)
