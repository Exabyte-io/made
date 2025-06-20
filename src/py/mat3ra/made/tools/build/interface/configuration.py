from typing import Optional, Union, List

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from mat3ra.made.material import Material
from mat3ra.made.tools.build.vacuum.configuration import VacuumConfiguration
from .. import BaseConfiguration, BaseConfigurationPydantic
from ..slab.configuration import (
    SlabConfiguration,
    SlabStrainedSupercellConfiguration,
    SlabStrainedSupercellWithGapConfiguration,
)

from ...analyze.other import get_chemical_formula


class InterfaceConfiguration(BaseConfigurationPydantic):
    # components and their modifiers added in the order they are stacked, from bottom to top
    stack_components: List[
        Union[
            SlabConfiguration,
            SlabStrainedSupercellConfiguration,
            SlabStrainedSupercellWithGapConfiguration,
            VacuumConfiguration,
        ]
    ]
    direction: AxisEnum = AxisEnum.z
    xy_shift: List[float] = [0.0, 0.0]  # in Angstroms

    @property
    def substrate_configuration(self) -> Union[SlabConfiguration, SlabStrainedSupercellConfiguration, SlabStrainedSupercellWithGapConfiguration]:
        return self.stack_components[0]

    @property
    def film_configuration(self) -> Union[SlabConfiguration, SlabStrainedSupercellConfiguration, SlabStrainedSupercellWithGapConfiguration]:
        return self.stack_components[1]

    def _get_crystal_from_slab_config(self, slab_config) -> Material:
        """Extract crystal material from any slab configuration type."""
        if hasattr(slab_config, 'atomic_layers'):
            # SlabConfiguration
            return slab_config.atomic_layers.crystal
        elif hasattr(slab_config, 'stack_components') and slab_config.stack_components:
            # SlabStrainedSupercellConfiguration and SlabStrainedSupercellWithGapConfiguration
            # The first stack component should be AtomicLayersUniqueRepeatedConfiguration
            atomic_layers = slab_config.stack_components[0]
            return atomic_layers.crystal
        else:
            raise ValueError(f"Cannot extract crystal from slab configuration: {type(slab_config)}")

    @property
    def vacuum_configuration(self) -> VacuumConfiguration:
        if len(self.stack_components) > 2:
            return self.stack_components[2]
        # Extract the crystal material from the film configuration
        film_crystal = self._get_crystal_from_slab_config(self.film_configuration)
        return VacuumConfiguration(size=0.0, crystal=film_crystal, direction=self.direction)

    @property
    def name(self) -> str:
        film_crystal = self._get_crystal_from_slab_config(self.film_configuration)
        substrate_crystal = self._get_crystal_from_slab_config(self.substrate_configuration)
        
        film_formula = get_chemical_formula(film_crystal)
        substrate_formula = get_chemical_formula(substrate_crystal)
        
        # Extract miller indices - handle both configuration types
        if hasattr(self.film_configuration, 'atomic_layers'):
            film_miller_indices = "".join([str(i) for i in self.film_configuration.atomic_layers.miller_indices])
        else:
            # For SlabStrainedSupercellConfiguration, get from first stack component
            film_miller_indices = "".join([str(i) for i in self.film_configuration.stack_components[0].miller_indices])
            
        if hasattr(self.substrate_configuration, 'atomic_layers'):
            substrate_miller_indices = "".join([str(i) for i in self.substrate_configuration.atomic_layers.miller_indices])
        else:
            # For SlabStrainedSupercellConfiguration, get from first stack component
            substrate_miller_indices = "".join([str(i) for i in self.substrate_configuration.stack_components[0].miller_indices])
            
        return f"{film_formula}({film_miller_indices})-{substrate_formula}({substrate_miller_indices}), Interface"


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
