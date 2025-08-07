from typing import Union

from mat3ra.made.material import Material
from .....build_components import MaterialWithBuildMetadata
from mat3ra.made.tools.build_components.entities.reusable.base_builder import BaseConfigurationPydantic


class MonolayerConfiguration(BaseConfigurationPydantic):
    """
    Configuration for creating monolayer structures.

    Miller indices are automatically determined based on crystal type:
    - HEX crystals: (0, 0, 1)
    - FCC/CUB crystals: (1, 1, 1) with primitive cell

    Args:
        crystal: The crystal material to create the monolayer from.
        vacuum: Size of the vacuum layer in Angstroms.
    """

    crystal: Union[Material, MaterialWithBuildMetadata]
    vacuum: float = 10.0
