from typing import List, Union, Any

from mat3ra.esse.models.materials_category_components.operations.core.combinations.stack import StackSchema
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from mat3ra.made.material import Material
from mat3ra.made.tools.build.vacuum.configuration import VacuumConfiguration


class StackConfiguration(StackSchema):

    stack_components: List[Union[Material, Any]]
    direction: AxisEnum = AxisEnum.z

    @property
    def atomic_layers(self):
        return self.stack_components[0]

    @property
    def vacuum_configuration(self) -> VacuumConfiguration:
        return self.stack_components[1]
