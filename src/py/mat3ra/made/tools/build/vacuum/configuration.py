from typing import Optional

from mat3ra.esse.models.materials_category_components.entities.core.two_dimensional.vacuum import (
    VacuumConfigurationSchema,
)

from mat3ra.made.material import Material
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from mat3ra.made.tools.build import BaseConfigurationPydantic


class VacuumConfiguration(VacuumConfigurationSchema, BaseConfigurationPydantic):
    type: str = "VacuumConfiguration"
    size: float = VacuumConfigurationSchema.model_fields["size"].default
    crystal: Material
    direction: Optional[AxisEnum] = VacuumConfigurationSchema.model_fields["direction"].default
