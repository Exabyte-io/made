from typing import Tuple

from pydantic import BaseModel

from mat3ra.made.material import Material
from .. import BaseConfiguration, BaseConfigurationPydantic


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
    
    crystal: Material
    vacuum: float = 10.0


