from typing import Optional

from mat3ra.code.entity import InMemoryEntity
from pydantic import BaseModel

from mat3ra.made.material import Material
from .termination_pair import TerminationPair
from .. import BaseConfiguration
from ..slab import Termination
from ..slab.configuration import SlabConfiguration


class InterfaceConfiguration(BaseModel, InMemoryEntity):
    """
    Configuration for an interface between two slabs.

    Args:
        film_configuration (SlabConfiguration): The configuration of the film slab.
        substrate_configuration (SlabConfiguration): The configuration of the substrate slab.
        film_termination (Termination): The termination of the film.
        substrate_termination (Termination): The termination of the substrate.
        distance_z (float): The distance between the film and the substrate in Angstroms.
        vacuum (float): The vacuum thickness, in Angstroms.
    """

    film_configuration: SlabConfiguration
    substrate_configuration: SlabConfiguration
    film_termination: Termination
    substrate_termination: Termination
    distance_z: float = 3.0
    vacuum: float = 5.0

    @property
    def termination_pair(self):
        return TerminationPair(self.film_termination, self.substrate_termination)

    @property
    def _json(self):
        return {
            "type": "InterfaceConfiguration",
            "film_configuration": self.film_configuration.to_json(),
            "substrate_configuration": self.substrate_configuration.to_json(),
            "film_termination": str(self.film_termination),
            "substrate_termination": str(self.substrate_termination),
            "distance_z": self.distance_z,
            "vacuum": self.vacuum,
        }


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
