from typing import List

import numpy as np
from pydantic import BaseModel, ConfigDict

from .bond_directions import BondDirections


class BondDirectionsTemplatesForElement(BaseModel):
    bond_directions_templates: List[BondDirections]
    element: str

    model_config = ConfigDict(arbitrary_types_allowed=True)

    def to_ndarray(self):
        return [np.array(template) for template in self.bond_directions_templates]
