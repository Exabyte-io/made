from typing import Optional, List, Union

from mat3ra.code.entity import InMemoryEntity
from pydantic import BaseModel

from mat3ra.made.material import Material
from mat3ra.made.utils import get_atomic_coordinates_extremum
from .enums import (
    PointDefectTypeEnum,
    AtomPlacementMethodEnum,
    ComplexDefectTypeEnum,
)
from ...analyze.other import get_closest_site_id_from_coordinate


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


class PointDefectConfigurationLegacy(BaseDefectConfiguration, InMemoryEntity):
    """
    Configuration for a point defect.

    Defect spatial arrangement can be set in the following ways:
    1. **Using Exact Coordinates**: Directly specify the exact coordinates of the defect in the crystal.
    2. **Using Approximate Coordinates**: Provide an approximate coordinate using from_approximate_coordinate method.
    3. **Using Site ID**: Specify the defect by providing a site ID using from_site_id method.


    Args:
        crystal (Material): The Material object.
        defect_type (PointDefectTypeEnum): The type of the defect.
        coordinate (Optional[List[float]]): The crystal coordinate of the defect.
        chemical_element (Optional[str]): The chemical element.
        use_cartesian_coordinates (bool): Whether to use cartesian coordinates for coordinates.
    """

    defect_type: PointDefectTypeEnum
    coordinate: List[float] = [0, 0, 0]
    chemical_element: Optional[str] = None
    use_cartesian_coordinates: bool = False
    placement_method: Optional[AtomPlacementMethodEnum] = None

    @classmethod
    def from_site_id(
        cls, crystal: Material, defect_type: PointDefectTypeEnum, site_id: int, chemical_element: Optional[str] = None
    ):
        if not crystal:
            raise RuntimeError("Crystal is not defined")
        coordinate = crystal.coordinates_array[site_id]
        return cls(crystal=crystal, defect_type=defect_type, coordinate=coordinate, chemical_element=chemical_element)

    @classmethod
    def from_approximate_coordinate(
        cls,
        crystal: Material,
        defect_type: PointDefectTypeEnum,
        approximate_coordinate: List[float],
        chemical_element: Optional[str] = None,
        use_cartesian_coordinates: bool = False,
    ):
        if not crystal:
            raise RuntimeError("Crystal is not defined")
        closest_site_id = get_closest_site_id_from_coordinate(
            material=crystal,
            coordinate=approximate_coordinate,
            use_cartesian_coordinates=use_cartesian_coordinates,
        )
        return cls.from_site_id(
            crystal=crystal, defect_type=defect_type, site_id=closest_site_id, chemical_element=chemical_element
        )

    @classmethod
    def from_dict(cls, crystal: Material, data: dict):
        if "site_id" in data:
            return cls.from_site_id(crystal=crystal, **data)
        elif "coordinate" in data:
            return cls(crystal=crystal, **data)
        elif "approximate_coordinate" in data:
            return cls.from_approximate_coordinate(crystal=crystal, **data)
        else:
            raise ValueError(f"Invalid defect configuration: {data}")

    @property
    def _json(self):
        return {
            **super()._json,
            "type": "PointDefectConfiguration",
            "defect_type": self.defect_type.name,
            "coordinate": self.coordinate,
            "chemical_element": self.chemical_element,
            "use_cartesian_coordinates": self.use_cartesian_coordinates,
            "placement_method": self.placement_method.name if self.placement_method else None,
        }


class SlabDefectConfigurationLegacy(BaseDefectConfiguration, InMemoryEntity):
    """
    Configuration for a slab defect.

    Args:
        crystal (Material): The Material object.
        number_of_added_layers (Union[int, float]): The number of added layers to the slab.
    """

    number_of_added_layers: Union[int, float] = 1

    @property
    def _json(self):
        return {
            **super()._json,
            "type": self.get_cls_name(),
            "number_of_added_layers": self.number_of_added_layers,
        }


class SlabPointDefectConfiguration(SlabDefectConfigurationLegacy, PointDefectConfigurationLegacy):
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
    placement_method: AtomPlacementMethodEnum = AtomPlacementMethodEnum.EXACT_COORDINATE

    @property
    def _json(self):
        return {
            **super()._json,
            "type": self.get_cls_name(),
            "defect_type": self.defect_type.name,
            "placement_method": self.placement_method.name,
        }


class PointDefectPairConfiguration(BaseDefectConfiguration, InMemoryEntity):
    """
    Configuration for a pair of point defects.

    Args:
        primary_defect_configuration: The first defect configuration.
        secondary_defect_configuration: The second defect configuration. Material is used from the primary defect.
    """

    defect_type: ComplexDefectTypeEnum = ComplexDefectTypeEnum.PAIR
    primary_defect_configuration: Union[PointDefectConfigurationLegacy, AdatomSlabPointDefectConfiguration]
    secondary_defect_configuration: Union[PointDefectConfigurationLegacy, AdatomSlabPointDefectConfiguration]

    @property
    def _json(self):
        return {
            "type": self.get_cls_name(),
            "defect_type": self.defect_type.name,
            "primary_defect_configuration": self.primary_defect_configuration.to_json(),
            "secondary_defect_configuration": self.secondary_defect_configuration.to_json(),
        }
