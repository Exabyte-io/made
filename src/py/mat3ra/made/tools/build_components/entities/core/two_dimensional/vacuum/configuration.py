from typing import Optional, Union

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.esse.models.materials_category_components.entities.core.two_dimensional.vacuum import (
    VacuumConfigurationSchema,
)
from mat3ra.made.material import Material
from mat3ra.made.tools.build_components import MaterialWithBuildMetadata
from mat3ra.made.tools.build_components.entities.reusable.three_dimensional.crystal_lattice_base.base_configuration_pydantic import (
    BaseConfigurationPydantic,
)


class VacuumConfiguration(VacuumConfigurationSchema, BaseConfigurationPydantic):
    type: str = "VacuumConfiguration"
    size: float = VacuumConfigurationSchema.model_fields["size"].default
    crystal: Union[Material, MaterialWithBuildMetadata, None] = None
    direction: Optional[AxisEnum] = AxisEnum.z
