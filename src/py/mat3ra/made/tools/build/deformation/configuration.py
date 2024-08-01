from typing import Callable, Dict, Tuple

from mat3ra.code.entity import InMemoryEntity
from mat3ra.made.material import Material
from pydantic import BaseModel

from ...utils import DeformationFunctionHolder


class DeformationConfiguration(BaseModel, InMemoryEntity):
    slab: Material
    deformation_function: Tuple[Callable, Dict] = DeformationFunctionHolder.sine_wave()
    use_cartesian_coordinates: bool = True

    class Config:
        arbitrary_types_allowed = True

    @property
    def _json(self):
        deformation_function_json = self.deformation_function[1]
        return {
            "type": self.get_cls_name(),
            "slab": self.slab.to_json(),
            "deformation_function": deformation_function_json,
            "use_cartesian_coordinates": self.use_cartesian_coordinates,
        }
