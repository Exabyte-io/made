from typing import Union

from mat3ra.made.material import Material
from .. import BaseConfigurationPydantic, MaterialWithBuildMetadata
from .functions import FunctionHolder, SineWavePerturbationFunctionHolder, PerturbationFunctionHolder


class PerturbationConfiguration(BaseConfigurationPydantic):
    """
    Configuration for a geometrical perturbation.

    Args:
        material (Material): The Material object.
        perturbation_function_holder (PerturbationFunctionHolder): The perturbation function holder.
        use_cartesian_coordinates (bool): Whether to use cartesian coordinates
    """

    material: Union[Material, MaterialWithBuildMetadata]
    perturbation_function_holder: Union[
        SineWavePerturbationFunctionHolder, PerturbationFunctionHolder, FunctionHolder
    ] = SineWavePerturbationFunctionHolder()
    use_cartesian_coordinates: bool = True
