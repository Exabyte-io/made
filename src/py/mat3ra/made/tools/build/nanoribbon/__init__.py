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
    # create a supercell of the material (must be 2D) with the desired width and length
    n = config.length // config.width + 1
    m = config.width // config.length + 1
    scaling_matrix = np.diag([config.width * n, config.length * m, 1])
    supercell = create_supercell(config.material, scaling_matrix)

    a_nudge = 0.1
    min_coordinate = [a_nudge, a_nudge, 0]
    max_coordinate = [
        config.width * config.material.lattice.a - a_nudge,
        config.length * config.material.lattice.b - a_nudge,
        config.material.lattice.c,
    ]

    nanoribbon = filter_by_rectangle_projection(
        supercell, min_coordinate=min_coordinate, max_coordinate=max_coordinate, use_cartesian_coordinates=True
    )

    return nanoribbon
