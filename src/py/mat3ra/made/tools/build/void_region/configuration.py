from typing import Union

from mat3ra.esse.models.materials_category_components.entities.auxiliary.zero_dimensional.void_region import (
    VoidRegionSchema,
)

from .. import MaterialWithBuildMetadata
from ...utils.coordinate import SphereCoordinateCondition, BoxCoordinateCondition, CoordinateCondition


class VoidRegionConfiguration(VoidRegionSchema):
    crystal: MaterialWithBuildMetadata
    coordinate_condition: Union[SphereCoordinateCondition, BoxCoordinateCondition, CoordinateCondition]
    use_cartesian_coordinates: bool = False
    invert_selection: bool = False

    @property
    def condition_name(self):
        return self.coordinate_condition.__class__.__name__.replace("CoordinateCondition", "")
