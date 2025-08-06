from typing import List

import numpy as np
from mat3ra.made.tools.bond_directions.bond_directions import BondDirections
from pydantic import BaseModel, ConfigDict


class BondDirectionsTemplatesForElement(BaseModel):
    bond_directions_templates: List[BondDirections]
    element: str

    model_config = ConfigDict(arbitrary_types_allowed=True)

    def to_ndarray(self):
        return [np.array(template) for template in self.bond_directions_templates]
