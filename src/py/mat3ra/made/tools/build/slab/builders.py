from mat3ra.made.material import Material
from .configuration import (
    SlabConfiguration,
)
from ..stack.builders import StackBuilder2Components
from ..vacuum.builders import VacuumBuilder
from ...build import BaseBuilder
from ...operations.core.unary import stack, supercell


class SlabBuilderParameters(BaseBuilder):
    make_primitive: bool = False
    use_orthogonal_c: bool = True


class SlabBuilder(StackBuilder2Components):

    def generate(self, configuration: SlabConfiguration) -> Material:
        stacked_materials = super().generate(configuration)
        supercell_slab = supercell(stacked_materials, configuration.xy_supercell_matrix)
        return supercell_slab

    def get_material(self, configuration: SlabConfiguration) -> Material:
        return self.generate(configuration)
