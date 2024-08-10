from typing import Literal

import numpy as np
from mat3ra.made.tools.build.supercell import create_supercell
from mat3ra.made.tools.modify import filter_by_condition_on_coordinates, rotate_material, filter_by_rectangle_projection
from pydantic import BaseModel

from mat3ra.made.material import Material


from pydantic import BaseModel


class NanoribbonConfiguration(BaseModel):
    material: Material
    width: int  # in number of unit cells
    length: int  # in number of unit cells
    edge_type: Literal["armchair", "zigzag"]

    class Config:
        arbitrary_types_allowed = True


def build_nanoribbon(config: NanoribbonConfiguration) -> Material:
    n = max(config.width, config.length)
    scaling_matrix = np.diag([n, n, 1])
    supercell = create_supercell(config.material, scaling_matrix)
    width_crystal = config.width / config.length
    min_coordinate = [0, 0, 0]
    max_coordinate = [1, width_crystal, 1]
    if config.edge_type == "armchair":
        rotation_matrix = [[1, -1, 0], [1, 1, 0], [0, 0, 1]]
        supercell = create_supercell(supercell, supercell_matrix=rotation_matrix)

    nanoribbon = filter_by_rectangle_projection(
        supercell, min_coordinate=min_coordinate, max_coordinate=max_coordinate, use_cartesian_coordinates=False
    )

    return nanoribbon
