from typing import Union

from mat3ra.made.material import Material
from .....build_components import MaterialWithBuildMetadata
from .builder import MonolayerBuilder
from .configurations import MonolayerConfiguration


def create_monolayer(
    crystal: Union[Material, MaterialWithBuildMetadata],
    vacuum: float = 10.0,
) -> Material:
    """
    Creates a monolayer material from a crystal material.

    Miller indices are automatically determined based on crystal type:
    - HEX crystals: (0, 0, 1) - basal plane
    - FCC/CUB crystals: (1, 1, 1) with primitive cell
    - Other types: (1, 1, 1) as generic fallback

    The function creates different monolayer structures based on the crystal type:
    - HEX: Creates a slab with thickness=1, translates to center, then filters half
    - FCC/CUB: Creates a slab with thickness=1 using primitive cell, then applies
      specific filtering and centering operations

    Args:
        crystal: The crystal material to create the monolayer from.
        vacuum: Size of the vacuum layer in Angstroms (default: 10.0).

    Returns:
        Material: The generated monolayer material.

    """
    configuration = MonolayerConfiguration(
        crystal=crystal,
        vacuum=vacuum,
    )

    builder = MonolayerBuilder()
    return builder.get_material(configuration)
