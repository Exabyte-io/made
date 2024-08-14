from typing import Literal, List, Optional, Any

import numpy as np
from mat3ra.code.entity import InMemoryEntity
from mat3ra.made.tools.build import BaseBuilder
from mat3ra.made.tools.build.supercell import create_supercell
from mat3ra.made.tools.modify import filter_by_rectangle_projection, wrap_to_unit_cell
from pydantic import BaseModel

from mat3ra.made.material import Material


class NanoribbonConfiguration(BaseModel, InMemoryEntity):
    """
    Configuration for building a nanoribbon.

    Attributes:
        material (Material): The material to build the nanoribbon from.
        width (int): The width of the nanoribbon in number of unit cells.
        length (int): The length of the nanoribbon in number of unit cells.
        edge_type (Literal["armchair", "zigzag"]): The edge type of the nanoribbon.
    """

    material: Material
    width: int  # in number of unit cells
    length: int  # in number of unit cells
    vacuum_width: int = 3  # in number of unit cells
    vacuum_length: int = 0  # in number of unit cells
    edge_type: Literal["armchair", "zigzag"]

    class Config:
        arbitrary_types_allowed = True

    @property
    def _json(self):
        return {
            "material": self.material.to_json(),
            "width": self.width,
            "length": self.length,
            "vacuum_width": self.vacuum_width,
            "vacuum_length": self.vacuum_length,
            "edge_type": self.edge_type,
        }


class NanoribbonBuilder(BaseBuilder):
    _ConfigurationType: type(NanoribbonConfiguration) = NanoribbonConfiguration  # type: ignore
    _GeneratedItemType: Material = Material
    _PostProcessParametersType: Any = None

    def build_nanoribbon(self, config: NanoribbonConfiguration) -> Material:
        scaling_matrix = np.diag([config.length + config.vacuum_length, config.width + config.vacuum_width, 1])
        supercell = create_supercell(config.material, scaling_matrix)
        n = config.length + config.vacuum_length
        m = config.width + config.vacuum_width
        nudge_value = 0.001

        if config.edge_type == "armchair":
            rotation_matrix = [[1, -1, 0], [1, 1, 0], [0, 0, 1]]
            n = m = max(config.length + config.vacuum_length, config.width + config.vacuum_width)
            scaling_matrix = np.diag([n, n, 1])
            supercell = create_supercell(config.material, scaling_matrix)
            supercell = create_supercell(supercell, supercell_matrix=rotation_matrix)
        length_crystal = config.length / n
        width_crystal = config.width / m
        min_coordinate = [(1 - length_crystal) / 2, (1 - width_crystal) / 2 + nudge_value, 0]
        max_coordinate = [(1 + length_crystal) / 2 - nudge_value, (1 + width_crystal) / 2, 1]
        nanoribbon = filter_by_rectangle_projection(
            supercell, min_coordinate=min_coordinate, max_coordinate=max_coordinate, use_cartesian_coordinates=False
        )

        return nanoribbon

    def _generate(self, configuration: NanoribbonConfiguration) -> List[_GeneratedItemType]:
        nanoribbon = self.build_nanoribbon(configuration)
        return [nanoribbon]

    def _post_process(
        self,
        items: List[_GeneratedItemType],
        post_process_parameters: Optional[_PostProcessParametersType],
    ) -> List[Material]:
        return [wrap_to_unit_cell(item) for item in items]

    def _update_material_name(self, material: Material, configuration: NanoribbonConfiguration) -> Material:
        edge_type = configuration.edge_type.capitalize()
        material.name = f"{material.name} ({edge_type} nanoribbon)"
        return material
