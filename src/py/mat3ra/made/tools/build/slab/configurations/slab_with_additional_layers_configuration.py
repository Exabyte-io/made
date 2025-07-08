from typing import Union

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.esse.models.materials_category.pristine_structures.two_dimensional.slab import SlabConfigurationSchema

from mat3ra.made.material import Material
from .slab_configuration import SlabConfiguration
from ...stack.configuration import StackConfiguration
from ...vacuum.configuration import VacuumConfiguration


class SlabWithAdditionalLayersConfiguration(SlabConfigurationSchema, StackConfiguration):
    type: str = "SlabWithAdditionalLayersConfiguration"
    stack_components: list  # Same as SlabConfiguration
    direction: AxisEnum = AxisEnum.z
    number_of_additional_layers: Union[int, float] = 1

    @property
    def atomic_layers(self):
        return self.stack_components[0]

    @property
    def vacuum_configuration(self) -> VacuumConfiguration:
        return self.stack_components[1]

    @classmethod
    def from_parameters(
        cls,
        material_or_dict: Union[Material, dict],
        miller_indices: tuple,
        number_of_layers: int,
        number_of_additional_layers: Union[int, float] = 1,
        termination_formula: str = None,
        vacuum: float = 10.0,
    ) -> "SlabWithAdditionalLayersConfiguration":
        base_slab_configuration = SlabConfiguration.from_parameters(
            material_or_dict=material_or_dict,
            miller_indices=miller_indices,
            number_of_layers=number_of_layers,
            termination_formula=termination_formula,
            vacuum=vacuum,
        )
        return cls(
            stack_components=base_slab_configuration.stack_components,
            direction=base_slab_configuration.direction,
            number_of_additional_layers=number_of_additional_layers,
        )
