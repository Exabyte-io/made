from typing import Optional, List, Union, Generic, TypeVar

from mat3ra.code.entity import InMemoryEntity
from mat3ra.esse.models.materials_category.defects.one_dimensional.terrace.configuration import (
    TerraceConfigurationSchema,
)
from mat3ra.esse.models.materials_category.defects.slab.configuration import SlabDefectConfigurationSchema
from mat3ra.esse.models.materials_category.defects.zero_dimensional.adatom.configuration import (
    AdatomConfigurationSchema,
)
from mat3ra.esse.models.materials_category.defects.zero_dimensional.defect_pair.configuration import (
    PointDefectPairConfigurationSchema,
)
from mat3ra.esse.models.materials_category.defects.zero_dimensional.slab.configuration import (
    SlabPointDefectConfigurationSchema,
)
from mat3ra.esse.models.materials_category.single_material.defect.point_defect.configuration import (
    PointDefectConfigurationSchema,
)

from mat3ra.made.material import Material
from .enums import (
    PointDefectTypeEnum,
    SlabDefectTypeEnum,
    CoordinatesShapeEnum,
)
from .. import BaseConfigurationPydantic
from ...analyze.other import get_closest_site_id_from_coordinate, get_atomic_coordinates_extremum
from ...utils.coordinate import (
    CylinderCoordinateCondition,
    SphereCoordinateCondition,
    BoxCoordinateCondition,
    TriangularPrismCoordinateCondition,
    PlaneCoordinateCondition,
    CoordinateCondition,
)


class BaseDefectConfiguration(BaseConfigurationPydantic):
    type: str = "BaseDefectConfiguration"
    crystal: Material = None


class PointDefectConfiguration(PointDefectConfigurationSchema, BaseDefectConfiguration):
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

    type: str = "PointDefectConfiguration"
    crystal: Material = None

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


class SlabDefectConfiguration(BaseDefectConfiguration, InMemoryEntity):
    """
    Configuration for a slab defect.

    Args:
        crystal (Material): The Material object.
        number_of_added_layers (Union[int, float]): The number of added layers to the slab.
    """

    number_of_added_layers: Union[int, float] = 1


class SlabDefectConfigurationPydantic(SlabDefectConfigurationSchema, BaseDefectConfiguration):
    """
    Configuration for a slab defect.

    Args:
        crystal (Material): The Material object.
        number_of_added_layers (Union[int, float]): The number of added layers to the slab.
    """

    type: str = "SlabDefectConfiguration"
    crystal: Material = None


class SlabPointDefectConfiguration(
    SlabPointDefectConfigurationSchema, SlabDefectConfigurationPydantic, PointDefectConfiguration
):
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

    type: str = "SlabPointDefectConfiguration"

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


class AdatomSlabPointDefectConfiguration(AdatomConfigurationSchema, SlabPointDefectConfiguration):
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

    type: str = "AdatomSlabPointDefectConfiguration"


CoordinateConditionType = TypeVar("CoordinateConditionType", bound=CoordinateCondition)


class IslandSlabDefectConfiguration(SlabDefectConfigurationPydantic, Generic[CoordinateConditionType]):
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
        CoordinateConditionType,
    ] = CylinderCoordinateCondition()

    @classmethod
    def from_dict(cls, crystal: Material, condition: dict, **kwargs):
        """
        Creates an IslandSlabDefectConfiguration instance from a dictionary.

        Args:
            crystal (Material): The material object.
            condition (dict): The dictionary with shape and other parameters for the condition.
            kwargs: Other configuration parameters (like number_of_added_layers).
        """
        condition_obj = cls.get_coordinate_condition(shape=condition["shape"], dict_params=condition)
        return cls(crystal=crystal, condition=condition_obj, **kwargs)

    @staticmethod
    def get_coordinate_condition(shape: CoordinatesShapeEnum, dict_params: dict):
        """
        Returns the appropriate coordinate condition based on the shape provided.

        Args:
            shape (CoordinatesShapeEnum): Shape of the island (e.g., cylinder, box, etc.).
            dict_params (dict): Parameters for the shape condition.

        Returns:
            CoordinateCondition: The appropriate condition object.
        """
        if shape == "cylinder":
            return CylinderCoordinateCondition(
                center_position=dict_params.get("center_position", [0.5, 0.5]),
                radius=dict_params["radius"],
                min_z=dict_params["min_z"],
                max_z=dict_params["max_z"],
            )
        elif shape == "sphere":
            return SphereCoordinateCondition(
                center_position=dict_params.get("center_position", [0.5, 0.5]), radius=dict_params["radius"]
            )
        elif shape == "box":
            return BoxCoordinateCondition(
                min_coordinate=dict_params["min_coordinate"], max_coordinate=dict_params["max_coordinate"]
            )
        elif shape == "triangular_prism":
            return TriangularPrismCoordinateCondition(
                position_on_surface_1=dict_params["position_on_surface_1"],
                position_on_surface_2=dict_params["position_on_surface_2"],
                position_on_surface_3=dict_params["position_on_surface_3"],
                min_z=dict_params["min_z"],
                max_z=dict_params["max_z"],
            )
        else:
            raise ValueError(f"Unsupported island shape: {shape}")

    @property
    def _json(self):
        return {
            **super()._json,
            "type": self.get_cls_name(),
            "defect_type": self.defect_type.name,
            "condition": self.condition.get_json(),
        }


class TerraceSlabDefectConfiguration(TerraceConfigurationSchema, SlabDefectConfigurationPydantic):
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


class PointDefectPairConfiguration(PointDefectPairConfigurationSchema, BaseDefectConfiguration):
    """
    Configuration for a pair of point defects.

    Args:
        primary_defect_configuration: The first defect configuration.
        secondary_defect_configuration: The second defect configuration. Material is used from the primary defect.
    """

    type: str = "PointDefectPairConfiguration"
