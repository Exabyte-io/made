from typing import Optional, List, Any, Callable
from pydantic import BaseModel

from mat3ra.code.entity import InMemoryEntity
from mat3ra.made.material import Material

from ...analyze import get_closest_site_id_from_coordinate, get_atomic_coordinates_extremum
from .enums import PointDefectTypeEnum, SlabDefectTypeEnum
from ...utils import is_coordinate_in_cylinder


class BaseDefectConfiguration(BaseModel):
    # TODO: fix arbitrary_types_allowed error and set Material class type
    crystal: Any = None


class PointDefectConfiguration(BaseDefectConfiguration, InMemoryEntity):
    defect_type: PointDefectTypeEnum
    coordinate: List[float] = [0, 0, 0]  # fractional coordinates
    chemical_element: Optional[str] = None

    @classmethod
    def from_site_id(
        cls, crystal: Material, defect_type: PointDefectTypeEnum, site_id: int, chemical_element: Optional[str] = None
    ):
        if not crystal:
            raise RuntimeError("Crystal is not defined")
        coordinate = crystal.coordinates_array[site_id]
        return cls(crystal=crystal, defect_type=defect_type, coordinate=coordinate, chemical_element=chemical_element)

    @classmethod
    def from_approximate_position(
        cls,
        crystal: Material,
        defect_type: PointDefectTypeEnum,
        approximate_coordinate: List[float],
        chemical_element: Optional[str] = None,
    ):
        if not crystal:
            raise RuntimeError("Crystal is not defined")
        closest_site_id = get_closest_site_id_from_coordinate(crystal, approximate_coordinate)
        return cls.from_site_id(
            crystal=crystal, defect_type=defect_type, site_id=closest_site_id, chemical_element=chemical_element
        )

    @property
    def _json(self):
        return {
            "type": "PointDefectConfiguration",
            "crystal": self.crystal.to_json(),
            "defect_type": self.defect_type.name,
            "coordinate": self.coordinate,
            "chemical_element": self.chemical_element,
        }


class SlabDefectConfiguration(BaseDefectConfiguration, InMemoryEntity):
    pass


class SlabPointDefectConfiguration(SlabDefectConfiguration, PointDefectConfiguration):
    position_on_surface: List[float]
    distance_z: float

    def __init__(self, **data):
        super().__init__(**data)
        if not self.position_on_surface:
            self.position_on_surface = [self.coordinate[0], self.coordinate[1]]
        if not self.distance_z:
            self.distance_z = self.crystal.basis.cell.convert_point_to_cartesian(
                [
                    self.coordinate[0],
                    self.coordinate[1],
                    self.coordinate[2] - get_atomic_coordinates_extremum(self.crystal),
                ]
            )[2]

    @property
    def _json(self):
        return {
            **super()._json,
            "type": "SlabPointDefectConfiguration",
            "position_on_surface": self.position_on_surface,
            "distance_z": self.distance_z,
        }


class AdatomSlabPointDefectConfiguration(SlabPointDefectConfiguration):
    defect_type: PointDefectTypeEnum = PointDefectTypeEnum.ADATOM

    @property
    def _json(self):
        return {
            **super()._json,
            "type": "AdatomSlabPointDefectConfiguration",
        }


class IslandSlabDefectConfiguration(SlabDefectConfiguration):
    def _default_condition(coordinate: List[float]) -> bool:
        return is_coordinate_in_cylinder(coordinate, [0.5, 0.5], radius=0.25)

    defect_type: SlabDefectTypeEnum = SlabDefectTypeEnum.ISLAND
    condition: Optional[Callable[[List[float]], bool]] = _default_condition
    thickness: int = 1

    @property
    def _json(self):
        return {
            "type": "IslandSlabDefectConfiguration",
            "crystal": self.crystal.to_json(),
            "defect_type": self.defect_type.name,
            "condition": self.condition,
            "thickness": self.thickness,
        }
