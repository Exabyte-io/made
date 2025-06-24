from typing import Optional, Union, List

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.esse.models.materials_category.compound_pristine_structures.two_dimensional.interface.configuration import (  # noqa: E501
    InterfaceConfigurationSchema,
)

from mat3ra.made.material import Material
from mat3ra.made.tools.build.vacuum.configuration import VacuumConfiguration
from .. import BaseConfiguration, BaseConfigurationPydantic
from ..slab.configurations import (
    SlabConfiguration,
    SlabStrainedSupercellConfiguration,
    SlabStrainedSupercellWithGapConfiguration,
)


class InterfaceConfiguration(InterfaceConfigurationSchema, BaseConfigurationPydantic):
    # components and their modifiers added in the order they are stacked, from bottom to top
    stack_components: List[
        Union[
            SlabStrainedSupercellConfiguration,
            SlabStrainedSupercellWithGapConfiguration,
            VacuumConfiguration,
        ]
    ]
    direction: AxisEnum = AxisEnum.z
    xy_shift: List[float] = InterfaceConfigurationSchema.model_fields["xy_shift"].default  # in Angstroms

    @property
    def substrate_configuration(self) -> SlabConfiguration:
        return self.stack_components[0]

    @property
    def film_configuration(self) -> SlabConfiguration:
        return self.stack_components[1]

    @property
    def vacuum_configuration(self) -> VacuumConfiguration:
        if len(self.stack_components) > 2:
            return self.stack_components[2]
        return VacuumConfiguration(
            size=0.0, crystal=self.film_configuration.atomic_layers.crystal, direction=self.direction
        )


class TwistedInterfaceConfiguration(BaseConfiguration):
    """
    Configuration for creating a twisted interface between two slabs with specified twist angle.

    Args:
        film (Material): The film material.
        substrate (Material): The substrate material.
        twist_angle (float): Twist angle in degrees.
        distance_z (float): Vertical distance between layers in Angstroms.
        vacuum (float): Vacuum thickness, in Angstroms.
    """

    film: Material
    substrate: Optional[Material] = None
    twist_angle: float = 0.0
    distance_z: float = 3.0
    vacuum: float = 0.0

    @property
    def _json(self):
        return {
            "type": self.get_cls_name(),
            "film": self.film.to_json(),
            "substrate": self.substrate.to_json() if self.substrate else None,
            "twist_angle": self.twist_angle,
            "distance_z": self.distance_z,
        }


class NanoRibbonTwistedInterfaceConfiguration(TwistedInterfaceConfiguration):
    """
    Configuration for creating a twisted interface between two nano ribbons with specified twist angle.

    Args:
        film (Material): The film material.
        substrate (Material): The substrate material.
        twist_angle (float): Twist angle in degrees.
        ribbon_width (int): Width of the nanoribbon in unit cells.
        ribbon_length (int): Length of the nanoribbon in unit cells.
        distance_z (float): Vertical distance between layers in Angstroms.
        vacuum_x (float): Vacuum along x on both sides, in Angstroms.
        vacuum_y (float): Vacuum along y on both sides, in Angstroms.
    """

    ribbon_width: int = 1
    ribbon_length: int = 1
    vacuum_x: float = 5.0
    vacuum_y: float = 5.0

    @property
    def _json(self):
        return {
            **super()._json,
            "type": self.get_cls_name(),
            "ribbon_width": self.ribbon_width,
            "ribbon_length": self.ribbon_length,
            "vacuum_x": self.vacuum_x,
            "vacuum_y": self.vacuum_y,
        }
