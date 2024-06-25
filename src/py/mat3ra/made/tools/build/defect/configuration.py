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
    position: List[float] = [0, 0, 0]  # fractional coordinates
    chemical_element: Optional[str] = None

    @classmethod
    def from_site_id(cls, site_id: int, crystal: Material):
        if not crystal:
            RuntimeError("Crystal is not defined")
        position = crystal.coordinates_array[site_id]            
        return cls(crystal=crystal, position=position)

    @classmethod
    def from_approximate_position(
        cls,
        approximate_position: List[float],
        crystal: Material
    ):
        if not crystal:
            RuntimeError("Crystal is not defined")
        closest_site_id = get_closest_site_id_from_position(crystal, approximate_position)
        return cls.from_site_id(closest_site_id, crystal)

    @property
    def _json(self):
        return {
            "type": "PointDefectConfiguration",
            "defect_type": self.defect_type.name,
            "position": self.position,
            "chemical_element": self.chemical_element,
        }
