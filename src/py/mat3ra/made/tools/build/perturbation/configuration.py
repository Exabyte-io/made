from typing import Callable, Dict, Tuple

from mat3ra.code.entity import InMemoryEntity
from mat3ra.made.material import Material
from pydantic import BaseModel

from ...utils import PerturbationFunctionHolder


class PerturbationConfiguration(BaseModel, InMemoryEntity):
    material: Material
    perturbation_function: Tuple[Callable, Dict] = PerturbationFunctionHolder.sine_wave()
    use_cartesian_coordinates: bool = True

    class Config:
        arbitrary_types_allowed = True

    @property
    def _json(self):
        perturbation_function_json = self.perturbation_function[1]
        return {
            "type": self.get_cls_name(),
            "material": self.material.to_json(),
            "perturbation_function": perturbation_function_json,
            "use_cartesian_coordinates": self.use_cartesian_coordinates,
        }
