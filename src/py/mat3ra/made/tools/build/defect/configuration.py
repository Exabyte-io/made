from typing import Optional, List, Union, Any

import numpy as np
from pydantic import BaseModel

from .enums import PointDefectTypeEnum


class BaseDefectConfiguration(BaseModel):
    # TODO: fix arbitrary_types_allowed error and set Material class type
    crystal: Any = None
    # defect_type type can be an Enum for a specific defect class (for point defect, 2d defect, etc.)
    defect_type: Union[PointDefectTypeEnum, None] = None


class BasePointDefectConfiguration(BaseDefectConfiguration):
    position: Optional[List[float]] = None  # fractional coordinates
    site_id: Optional[int] = 0

    def __init__(self, **data):
        super().__init__(**data)
        if self.position:
            self.site_id = self._get_site_id_from_position()

    # TODO: move to analysis module?
    def _get_site_id_from_position(self):
        closest_site_id = self.site_id
        if self.crystal is not None:
            coordinates = np.array(self.crystal.basis["coordinates"])
            position = np.array(self.position)
            distances = np.linalg.norm(coordinates - position, axis=1)
            closest_site_id = np.argmin(distances)
        return closest_site_id


class VacancyConfiguration(BasePointDefectConfiguration):
    defect_type: PointDefectTypeEnum = PointDefectTypeEnum.VACANCY


class SubstitutionConfiguration(BasePointDefectConfiguration):
    defect_type: PointDefectTypeEnum = PointDefectTypeEnum.SUBSTITUTION
    element: Optional[str] = "Si"


class InterstitialConfiguration(BasePointDefectConfiguration):
    defect_type: PointDefectTypeEnum = PointDefectTypeEnum.INTERSTITIAL
    element: Optional[str] = "Si"
