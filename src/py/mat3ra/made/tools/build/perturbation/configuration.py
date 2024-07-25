from typing import Callable

from mat3ra.code.entity import InMemoryEntity
from mat3ra.made.material import Material
from pydantic import BaseModel


class PerturbationConfiguration(BaseModel, InMemoryEntity):
    slab: Material
    perturbation_func: Callable[[float, float], float]
    amplitude: float

    class Config:
        arbitrary_types_allowed = True

    @property
    def _json(self):
        return {
            "type": "PerturbationConfiguration",
            "slab": self.slab.to_json(),
            "amplitude": self.amplitude,
        }
