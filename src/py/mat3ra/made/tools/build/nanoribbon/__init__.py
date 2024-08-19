from typing import Literal, List, Optional, Any

import numpy as np
from mat3ra.code.entity import InMemoryEntity
from mat3ra.made.lattice import Lattice
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
        material = config.material
        provided_width = config.width
        provided_length = config.length
        vacuum_w = config.vacuum_width
        vacuum_l = config.vacuum_length
        edge_type = config.edge_type
        coeff = 1
        if edge_type == "armchair":
            provided_length, provided_width = provided_width, provided_length
            vacuum_w, vacuum_l = vacuum_l, vacuum_w
            coeff = coeff * (provided_width // provided_length + 1)
        supercell = create_supercell(material, np.diag([2 * provided_length * coeff, 2 * provided_width, 1]))
        visualize(supercell, repetitions=[1, 1, 1], rotation="0x")
        nudge_value = 0.01
        conditional_nudge_value = nudge_value * (-1 * (edge_type == "armchair") + 1 * (edge_type == "zigzag"))

        length = provided_length * np.dot(np.array(material.basis.cell.vector1), np.array([1, 0, 0]))
        width = provided_width * np.dot(np.array(material.basis.cell.vector2), np.array([0, 1, 0]))
        height = np.dot(np.array(material.basis.cell.vector3), np.array([0, 0, 1]))

        vacuum_length = vacuum_l * np.dot(np.array(material.basis.cell.vector1), np.array([1, 0, 0]))
        vacuum_width = vacuum_w * np.dot(np.array(material.basis.cell.vector2), np.array([0, 1, 0]))

        min_coordinate = [-nudge_value, conditional_nudge_value, 0]
        max_coordinate = [length - nudge_value, width + conditional_nudge_value, height]

        nanoribbon = filter_by_rectangle_projection(
            supercell, min_coordinate=min_coordinate, max_coordinate=max_coordinate, use_cartesian_coordinates=True
        )
        new_basis = nanoribbon.basis.copy()
        new_basis.to_cartesian()
        new_length = length + vacuum_length
        new_width = width + vacuum_width

        length_vector = [new_length, 0, 0]
        width_vector = [0, new_width, 0]
        height_vector = [0, 0, height]
        if edge_type == "armchair":
            length_vector, width_vector = width_vector, length_vector

        new_basis.cell.vector1 = length_vector
        new_basis.cell.vector2 = width_vector
        new_basis.cell.vector3 = height_vector
        new_basis.to_crystal()
        nanoribbon.basis = new_basis

        lattice = Lattice.from_vectors_array([length_vector, width_vector, height_vector])
        nanoribbon.lattice = lattice

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
