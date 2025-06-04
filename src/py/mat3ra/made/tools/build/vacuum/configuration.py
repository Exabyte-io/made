from typing import Optional

from mat3ra.esse.models.materials_category_components.entities.core.two_dimensional.vacuum import (
    VacuumConfigurationSchema,
)

from mat3ra.made.material import Material
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum


class VacuumConfiguration(VacuumConfigurationSchema):
    size: float
    crystal: Material
    direction: Optional[AxisEnum]
