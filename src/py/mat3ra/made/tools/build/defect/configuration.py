from typing import Optional, List, Any
from pydantic import BaseModel

from mat3ra.code.entity import InMemoryEntity
from mat3ra.made.material import Material

from ...analyze import get_closest_site_id_from_position
from .enums import PointDefectTypeEnum


class BaseDefectConfiguration(BaseModel):
    # TODO: fix arbitrary_types_allowed error and set Material class type
    crystal: Any = None


class PointDefectConfiguration(BaseDefectConfiguration, InMemoryEntity):
    defect_type: PointDefectTypeEnum
    position: Optional[List[float]] = [0, 0, 0]  # fractional coordinates
    site_id: Optional[int] = None
    chemical_element: Optional[str] = None

    def __init__(self, position=position, site_id=None, **data):
        super().__init__(**data)
        if site_id is not None:
            self.position = self.crystal.coordinates_array[site_id]
        else:
            self.site_id = get_closest_site_id_from_position(self.crystal, position)

    @classmethod
    def from_site_id(cls, site_id: int, **data):
        crystal: Material = data.get("crystal")
        position = crystal.coordinates_array[site_id]
        return cls(position=position, **data)

    @property
    def _json(self):
        return {
            "defect_type": self.defect_type,
            "position": self.position,
            "chemical_element": self.chemical_element,
        }
