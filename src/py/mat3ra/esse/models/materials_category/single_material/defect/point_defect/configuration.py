from typing import List, Optional

from mat3ra.esse.models.materials_category.single_material.defect.base import BaseDefectConfigurationSchema
from mat3ra.esse.models.materials_category.single_material.defect.enums import (
    AtomPlacementMethodEnum,
    PointDefectTypeEnum,
)


class PointDefectConfigurationSchema(BaseDefectConfigurationSchema):
    """
    Schema for point defect configuration.

    Args:
        defect_type (PointDefectTypeEnum): The type of the defect.
        coordinate (List[float]): The crystal coordinate of the defect.
        chemical_element (Optional[str]): The chemical element.
        use_cartesian_coordinates (bool): Whether to use cartesian coordinates for coordinates.
        placement_method (Optional[AtomPlacementMethodEnum]): The method used to place the atom.
    """

    defect_type: PointDefectTypeEnum
    coordinate: List[float] = [0, 0, 0]
    chemical_element: Optional[str] = None
    use_cartesian_coordinates: bool = False
    placement_method: Optional[AtomPlacementMethodEnum] = None
