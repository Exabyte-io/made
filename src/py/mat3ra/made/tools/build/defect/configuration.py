from typing import Optional, List, Union, Any

from pydantic import BaseModel

from .enums import PointDefectTypeEnum


class BaseDefectConfiguration(BaseModel):
    # TODO: fix arbitrary_types_allowed error and set Material class type
    material: Any
    # defect_type type can be an Enum for a specific defect class (for point defect, 2d defect, etc.)
    defect_type: Union[PointDefectTypeEnum, None] = None


class PointDefectConfiguration(BaseDefectConfiguration):
    defect_type: PointDefectTypeEnum = PointDefectTypeEnum.VACANCY
    # TODO: should come from Enum of elements
    specie: Optional[str] = None
    # TODO: import coordinate type from ESSE
    # only for interstitials
    position_shift: Optional[List[float]] = [0, 0, 0]
