from typing import List, Union, Optional, Tuple

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from mat3ra.made.material import Material
from ..lattice_lines.configuration import EdgeTypes, get_miller_indices_from_edge_type
from ..nanotape.builders import NanoTapeBuilder
from ..nanotape.configuration import NanoTapeConfiguration
from ..slab.entities import Termination
from ..stack.configuration import StackConfiguration
from ..vacuum.configuration import VacuumConfiguration


class NanoribbonConfiguration(StackConfiguration):
    """
    Configuration for building a nanoribbon from a nanotape.
    Nanoribbon = [NanoTape, vacuum] stacked on X direction.

    Args:
        stack_components: List of configuration objects for nanoribbon components.
        direction: Direction along which to stack components.
    """

    type: str = "NanoribbonConfiguration"
    stack_components: List[Union[NanoTapeConfiguration, VacuumConfiguration]]
    direction: AxisEnum = AxisEnum.x

    @property
    def nanotape(self):
        """Get the nanotape configuration component."""
        return self.stack_components[0]

    @property
    def vacuum_configuration(self) -> VacuumConfiguration:
        """Get the vacuum configuration component."""
        return self.stack_components[1]

    @classmethod
    def from_parameters(
        cls,
        material: Material,
        miller_indices_2d: Optional[Tuple[int, int]] = None,
        edge_type: Optional[EdgeTypes] = EdgeTypes.zigzag,
        width: int = 2,
        length: int = 2,
        vacuum_width: float = 10.0,
        vacuum_length: float = 10.0,
        termination: Optional[Termination] = None,
    ) -> "NanoribbonConfiguration":
        """
        Create a NanoribbonConfiguration from parameters.

        Args:
            material: The monolayer material to create the nanoribbon from.
            miller_indices_2d: The (u,v) Miller indices for the nanoribbon direction.
            edge_type: Edge type string ("zigzag"/"armchair"). Optional if miller_indices_2d is provided.
            width: The width of the nanoribbon in number of unit cells.
            length: The length of the nanoribbon in number of unit cells.
            vacuum_width: The width of the vacuum region in Angstroms (cartesian).
            vacuum_length: The length of the vacuum region in Angstroms (cartesian).
            termination: The termination to use for the nanoribbon. If None, uses default termination.

        Returns:
            NanoribbonConfiguration: The nanoribbon configuration.
        """
        if miller_indices_2d is None and edge_type is None:
            raise ValueError("Either miller_indices_2d or edge_type must be provided")
        if miller_indices_2d is None and edge_type is not None:
            miller_indices_2d = get_miller_indices_from_edge_type(edge_type)

        nanotape_config = NanoTapeConfiguration.from_parameters(
            material=material,
            miller_indices_2d=miller_indices_2d,
            width=width,
            length=length,
            termination=termination,
            vacuum_width=vacuum_width,
        )

        nanotape_builder = NanoTapeBuilder()
        nanotape_material = nanotape_builder.get_material(nanotape_config)

        vacuum_config = VacuumConfiguration(
            size=vacuum_length,
            crystal=nanotape_material,
            direction=AxisEnum.x,
        )

        return cls(
            stack_components=[nanotape_config, vacuum_config],
            direction=AxisEnum.x,
        )
