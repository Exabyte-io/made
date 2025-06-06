from typing import List, Any

from mat3ra.code.entity import InMemoryEntityPydantic
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.esse.models.materials_category_components.operations.core.combinations.stack import StackSchema

from mat3ra.made.tools.build.vacuum.configuration import VacuumConfiguration


class StackConfiguration(StackSchema, InMemoryEntityPydantic):

    stack_components: List[Any]  # Configuration objects only, no Materials
    direction: AxisEnum = AxisEnum.z

    @property
    def vacuum_configuration(self) -> VacuumConfiguration:
        return self.stack_components[1]
