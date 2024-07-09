from typing import Optional, List, Any
from pydantic import BaseModel

from mat3ra.code.entity import InMemoryEntity
from mat3ra.made.material import Material

from ...analyze import get_closest_site_id_from_position
from .enums import PointDefectTypeEnum, SlabDefectTypeEnum


class BaseDefectConfiguration(BaseModel):
    # TODO: fix arbitrary_types_allowed error and set Material class type
    crystal: Any = None


class PointDefectConfiguration(BaseDefectConfiguration, InMemoryEntity):
    defect_type: PointDefectTypeEnum
    position: List[float] = [0, 0, 0]  # fractional coordinates
    chemical_element: Optional[str] = None

    @classmethod
    def from_site_id(
        cls, crystal: Material, defect_type: PointDefectTypeEnum, site_id: int, chemical_element: Optional[str] = None
    ):
        if not crystal:
            raise RuntimeError("Crystal is not defined")
        position = crystal.coordinates_array[site_id]
        return cls(crystal=crystal, defect_type=defect_type, position=position, chemical_element=chemical_element)

    @classmethod
    def from_approximate_position(
        cls,
        crystal: Material,
        defect_type: PointDefectTypeEnum,
        approximate_position: List[float],
        chemical_element: Optional[str] = None,
    ):
        if not crystal:
            raise RuntimeError("Crystal is not defined")
        closest_site_id = get_closest_site_id_from_position(crystal, approximate_position)
        return cls.from_site_id(
            crystal=crystal, defect_type=defect_type, site_id=closest_site_id, chemical_element=chemical_element
        )

    @property
    def _json(self):
        return {
            "type": "PointDefectConfiguration",
            "crystal": self.crystal.to_json(),
            "defect_type": self.defect_type.name,
            "position": self.position,
            "chemical_element": self.chemical_element,
        }


class SlabDefectConfiguration(BaseDefectConfiguration, InMemoryEntity):
    pass


class AdatomSlabDefectConfiguration(SlabDefectConfiguration):
    defect_type: SlabDefectTypeEnum = SlabDefectTypeEnum.ADATOM
    position_on_surface: List[float] = [0.5, 0.5]
    distance_z: float = 2.0
    chemical_element: str = "Si"

    @property
    def _json(self):
        return {
            "type": "AdatomSlabDefectConfiguration",
            "crystal": self.crystal.to_json(),
            "defect_type": self.defect_type.name,
            "position_on_surface": self.position_on_surface,
            "distance_z": self.distance_z,
            "chemical_element": self.chemical_element,
        }
