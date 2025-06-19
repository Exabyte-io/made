from typing import List, Union, Optional

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.esse.models.materials_category.pristine_structures.two_dimensional.slab import SlabConfigurationSchema
from mat3ra.esse.models.materials_category.pristine_structures.two_dimensional.slab_strained_supercell import (
    SlabStrainedSupercellConfigurationSchema,
)
from mat3ra.esse.models.materials_category.pristine_structures.two_dimensional.slab_strained_supercell_with_gap import (
    SlabStrainedSupercellWithGapConfigurationSchema,
)
from mat3ra.esse.models.materials_category_components.entities.reusable.two_dimensional.crystal_lattice_planes import (
    CrystalLatticePlanesSchema,
)

from mat3ra.made.material import Material
from mat3ra.made.tools.build.slab.entities import Termination
from .. import BaseConfigurationPydantic
from ..stack.configuration import StackConfiguration
from ..vacuum.configuration import VacuumConfiguration


class CrystalLatticePlanesConfiguration(CrystalLatticePlanesSchema, BaseConfigurationPydantic):
    crystal: Material

    @property
    def in_plane_vectors(self):
        # Two vectors in the plane of the Miller indices
        return self.crystal.lattice.vector_arrays[:2, :2]


class AtomicLayersUnique(CrystalLatticePlanesConfiguration):
    pass


class AtomicLayersUniqueRepeatedConfiguration(AtomicLayersUnique):
    termination_top: Termination
    number_of_repetitions: int = 1


class SlabConfiguration(SlabConfigurationSchema, StackConfiguration):
    type: str = "SlabConfiguration"
    stack_components: List[
        Union[AtomicLayersUnique, AtomicLayersUniqueRepeatedConfiguration, VacuumConfiguration]  # No Materials!
    ]
    direction: AxisEnum = AxisEnum.z

    @property
    def atomic_layers(self):
        return self.stack_components[0]

    @property
    def vacuum_configuration(self) -> VacuumConfiguration:
        return self.stack_components[1]


class SlabStrainedSupercellConfiguration(SlabStrainedSupercellConfigurationSchema, SlabConfiguration):
    pass


class SlabStrainedSupercellWithGapConfiguration(SlabStrainedSupercellWithGapConfigurationSchema, SlabConfiguration):
    gap: Optional[float] = None  # If provided, the film is shifted to have it as smallest distance to the substrate.
