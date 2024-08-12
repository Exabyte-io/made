from typing import Optional, List, Union
from pydantic import BaseModel

from mat3ra.code.entity import InMemoryEntity
from mat3ra.made.material import Material

from ...analyze import get_closest_site_id_from_coordinate, get_atomic_coordinates_extremum
from ...utils.coordinate import (
    CylinderCoordinateCondition,
    SphereCoordinateCondition,
    BoxCoordinateCondition,
    TriangularPrismCoordinateCondition,
    PlaneCoordinateCondition,
)
from .enums import PointDefectTypeEnum, SlabDefectTypeEnum, AtomPlacementMethodEnum, ComplexDefectTypeEnum


class BaseDefectConfiguration(BaseModel):
    crystal: Material = None

    class Config:
        arbitrary_types_allowed = True

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
    """
    Configuration for a slab defect.

    Args:
        crystal (Material): The Material object.
        number_of_added_layers (int): The number of added layers.
    """

    number_of_added_layers: int = 1

    @property
    def _json(self):
        return {
            **super()._json,
            "type": self.get_cls_name(),
            "number_of_added_layers": self.number_of_added_layers,
        }


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
            "type": self.get_cls_name(),
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
    placement_method: AtomPlacementMethodEnum = AtomPlacementMethodEnum.COORDINATE

    @property
    def _json(self):
        return {
            **super()._json,
            "type": self.get_cls_name(),
            "defect_type": self.defect_type.name,
            "placement_method": self.placement_method.name,
        }


class IslandSlabDefectConfiguration(SlabDefectConfiguration):
    """
    Configuration for an island slab defect.

    Args:
        crystal (Material): The Material object.
        defect_type (SlabDefectTypeEnum): The type of the defect.
        condition (Optional[Tuple[Callable[[List[float]], bool], Dict]]): The condition on coordinates
            to shape the island. Defaults to a cylinder.
        number_of_added_layers (int): The number of added layers to the slab which will form the island.
    """

    defect_type: SlabDefectTypeEnum = SlabDefectTypeEnum.ISLAND
    condition: Union[
        CylinderCoordinateCondition,
        SphereCoordinateCondition,
        BoxCoordinateCondition,
        TriangularPrismCoordinateCondition,
        PlaneCoordinateCondition,
    ] = CylinderCoordinateCondition()

    @property
    def _json(self):
        return {
            **super()._json,
            "type": self.get_cls_name(),
            "defect_type": self.defect_type.name,
            "condition": self.condition.get_json(),
        }


class TerraceSlabDefectConfiguration(SlabDefectConfiguration):
    """
    Configuration for a terrace slab defect.

    Args:
        crystal (Material): The Material object (must be a created slab).
        defect_type (SlabDefectTypeEnum): The type of the defect.
        cut_direction (List[int]): The direction of the cut as lattice vector, can be thought as a normal to the plane
            that cuts the slab with added number of layers.
        pivot_coordinate (List[float]): The pivot coordinate: the point in the unit cell
            where the normal of the cut plane passes through.
        number_of_added_layers (int): The number of added layers to the slab which will form the terrace.
        use_cartesian_coordinates (bool): The flag to use cartesian coordinates for coordinates and vectors.
        rotate_to_match_pbc (bool): The flag to rotate the slab with a terrace to match periodic boundary conditions.
    """

    defect_type: SlabDefectTypeEnum = SlabDefectTypeEnum.TERRACE
    cut_direction: List[int] = [1, 0, 0]
    pivot_coordinate: List[float] = [0.5, 0.5, 0.5]
    use_cartesian_coordinates: bool = False
    rotate_to_match_pbc: bool = True

    @property
    def _json(self):
        return {
            **super()._json,
            "type": self.get_cls_name(),
            "defect_type": self.defect_type.name,
            "cut_direction": self.cut_direction,
            "pivot_coordinate": self.pivot_coordinate,
            "use_cartesian_coordinates": self.use_cartesian_coordinates,
            "rotate_to_match_pbc": self.rotate_to_match_pbc,
        }


class PointDefectPairConfiguration(BaseDefectConfiguration, InMemoryEntity):
    """
    Configuration for a pair of point defects.

    Args:
        primary_defect_configuration: The first defect configuration.
        secondary_defect_configuration: The second defect configuration. Material is used from the primary defect.
    """

    defect_type: ComplexDefectTypeEnum = ComplexDefectTypeEnum.PAIR
    primary_defect_configuration: Union[PointDefectConfiguration, AdatomSlabPointDefectConfiguration]
    secondary_defect_configuration: Union[PointDefectConfiguration, AdatomSlabPointDefectConfiguration]

    @property
    def _json(self):
        return {
            "type": self.get_cls_name(),
            "defect_type": self.defect_type.name,
            "primary_defect_configuration": self.primary_defect_configuration.to_json(),
            "secondary_defect_configuration": self.secondary_defect_configuration.to_json(),
        }
