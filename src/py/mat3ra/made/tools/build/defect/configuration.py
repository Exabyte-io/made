from typing import Optional, List, Any, Callable, Dict, Tuple
from pydantic import BaseModel

from mat3ra.code.entity import InMemoryEntity
from mat3ra.made.material import Material

from ...analyze import get_closest_site_id_from_coordinate, get_atomic_coordinates_extremum
from ...utils import CoordinateConditionBuilder
from .enums import PointDefectTypeEnum, SlabDefectTypeEnum


class BaseDefectConfiguration(BaseModel):
    # TODO: fix arbitrary_types_allowed error and set Material class type
    crystal: Any = None

    @property
    def _json(self):
        return {
            "type": "BaseDefectConfiguration",
            "crystal": self.crystal.to_json(),
        }


class PointDefectConfiguration(BaseDefectConfiguration, InMemoryEntity):
    """
    Configuration for a point defect.

    Args:
        crystal (Material): The Material object.
        defect_type (PointDefectTypeEnum): The type of the defect.
        coordinate (List[float]): The crystal coordinate of the defect.
        chemical_element (Optional[str]): The chemical element.
    """

    defect_type: PointDefectTypeEnum
    coordinate: List[float] = [0, 0, 0]  # crystal coordinates
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
            **super()._json,
            "type": "PointDefectConfiguration",
            "defect_type": self.defect_type.name,
            "coordinate": self.coordinate,
            "chemical_element": self.chemical_element,
        }


class SlabDefectConfiguration(BaseDefectConfiguration, InMemoryEntity):
    pass


class SlabPointDefectConfiguration(SlabDefectConfiguration, PointDefectConfiguration):
    """
    Configuration for a slab point defect.

    Args:
        crystal (Material): The Material object.
        defect_type (PointDefectTypeEnum): The type of the defect.
        coordinate (List[float]): The crystal coordinate of the defect.
        chemical_element (Optional[str]): The chemical element.
        position_on_surface (List[float]): The position on the surface in 2D crystal coordinates.
        distance_z (float): The distance in z direction in angstroms.
    """

    position_on_surface: List[float]
    distance_z: float  # in angstroms

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
    """
    Configuration for an adatom slab point defect.

    Args:
        crystal (Material): The Material object.
        defect_type (PointDefectTypeEnum): The type of the defect.
        coordinate (List[float]): The crystal coordinate of the defect.
        chemical_element (Optional[str]): The chemical element.
        position_on_surface (List[float]): The position on the surface in 2D crystal coordinates.
        distance_z (float): The distance in z direction in angstroms.
    """

    defect_type: PointDefectTypeEnum = PointDefectTypeEnum.ADATOM

    @property
    def _json(self):
        return {
            **super()._json,
            "type": "AdatomSlabPointDefectConfiguration",
        }


class IslandSlabDefectConfiguration(SlabDefectConfiguration):
    """
    Configuration for an island slab defect.

    Args:
        crystal (Material): The Material object.
        defect_type (SlabDefectTypeEnum): The type of the defect.
        condition (Optional[Tuple[Callable[[List[float]], bool], Dict]]): The condition on coordinates
        to shape the island. Defaults to a cylinder.
        thickness (int): The thickness of the defect in atomic layers.
    """

    defect_type: SlabDefectTypeEnum = SlabDefectTypeEnum.ISLAND
    condition: Optional[Tuple[Callable[[List[float]], bool], Dict]] = CoordinateConditionBuilder().cylinder()
    thickness: int = 1  # in atomic layers

    @property
    def _json(self):
        _, condition_json = self.condition
        return {
            **super()._json,
            "type": "IslandSlabDefectConfiguration",
            "defect_type": self.defect_type.name,
            "condition": condition_json,
            "thickness": self.thickness,
        }
