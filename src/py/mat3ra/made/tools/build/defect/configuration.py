from typing import Optional, List, Any

from pydantic import BaseModel

from ...analyze import get_closest_site_id_from_position


class BaseDefectConfiguration(BaseModel):
    # TODO: fix arbitrary_types_allowed error and set Material class type
    crystal: Any = None


class PointDefectConfiguration(BaseDefectConfiguration):
    position: Optional[List[float]] = [0, 0, 0]  # fractional coordinates
    chemical_element: Optional[str] = None

    def __init__(self, position=position, site_id=None, **data):
        super().__init__(**data)
        if site_id is not None:
            self.position = self.crystal.coordinates_array[site_id]
        else:
            self.site_id = get_closest_site_id_from_position(self.crystal, position)

    @classmethod
    def from_site_id(cls, site_id: int, **data):
        position = cls.crystal.coordinates_array[site_id]
        return cls(position=position, **data)
