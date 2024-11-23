from typing import Union

from mat3ra.code.entity import InMemoryEntity
from mat3ra.made.material import Material
from pydantic import BaseModel

from ...utils.perturbation import SineWavePerturbationFunctionHolder, PerturbationFunctionHolder


class PerturbationConfiguration(BaseModel, InMemoryEntity):
    """
    Configuration for a geometrical perturbation.

    Args:
        material (Material): The Material object.
        perturbation_function_holder (PerturbationFunctionHolder): The perturbation function holder.
        use_cartesian_coordinates (bool): Whether to use cartesian coordinates
    """

    material: Material
    perturbation_function_holder: Union[
        SineWavePerturbationFunctionHolder, PerturbationFunctionHolder
    ] = SineWavePerturbationFunctionHolder()
    use_cartesian_coordinates: bool = True

    class Config:
        arbitrary_types_allowed = True

    @property
    def _json(self):
        perturbation_function_json = self.perturbation_function_holder.get_json()
        return {
            "type": self.get_cls_name(),
            "material": self.material.to_json(),
            "perturbation_function": perturbation_function_json,
            "use_cartesian_coordinates": self.use_cartesian_coordinates,
        }
