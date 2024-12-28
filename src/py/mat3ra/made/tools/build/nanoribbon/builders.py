from typing import List, Optional, Any, Tuple

import numpy as np

from mat3ra.made.material import Material
from mat3ra.made.tools.build import BaseBuilder
from mat3ra.made.tools.build.supercell import create_supercell
from mat3ra.made.tools.modify import filter_by_rectangle_projection, wrap_to_unit_cell

from ...modify import translate_to_center, rotate
from .configuration import NanoribbonConfiguration
from .enums import EdgeTypes


class NanoribbonBuilder(BaseBuilder):
    """
    Builder class for creating a nanoribbon from a material.

    The process creates a supercell with large enough dimensions to contain the nanoribbon and then
    filters the supercell to only include the nanoribbon. The supercell is then centered and returned as the nanoribbon.
    The nanoribbon can have either Armchair or Zigzag edge. The edge is defined along the vector 1 of the material cell,
    which corresponds to [1,0,0] direction.
    """

    _ConfigurationType: type(NanoribbonConfiguration) = NanoribbonConfiguration  # type: ignore
    _GeneratedItemType: Material = Material
    _PostProcessParametersType: Any = None

    def create_nanoribbon(self, config: NanoribbonConfiguration) -> Material:
        material = config.material
        (
            length_cartesian,
            width_cartesian,
            height_cartesian,
            vacuum_length_cartesian,
            vacuum_width_cartesian,
        ) = self._calculate_cartesian_dimensions(config, material)
        length_lattice_vector, width_lattice_vector, height_lattice_vector = self._get_new_lattice_vectors(
            length_cartesian,
            width_cartesian,
            height_cartesian,
            vacuum_length_cartesian,
            vacuum_width_cartesian,
            config.edge_type,
        )
        n = max(config.length, config.width)
        large_supercell_to_cut = create_supercell(material, np.diag([2 * n, 2 * n, 1]))

        min_coordinate, max_coordinate = self._calculate_coordinates_of_cut(
            length_cartesian, width_cartesian, height_cartesian, config.edge_type
        )
        nanoribbon = filter_by_rectangle_projection(
            large_supercell_to_cut,
            min_coordinate=min_coordinate,
            max_coordinate=max_coordinate,
            use_cartesian_coordinates=True,
        )
        nanoribbon.set_new_lattice_vectors(length_lattice_vector, width_lattice_vector, height_lattice_vector)
        return translate_to_center(nanoribbon)

    @staticmethod
    def _calculate_cartesian_dimensions(config: NanoribbonConfiguration, material: Material):
        """
        Calculate the dimensions of the nanoribbon in the cartesian coordinate system.
        """
        nanoribbon_width = config.width
        nanoribbon_length = config.length
        vacuum_width = config.vacuum_width
        vacuum_length = config.vacuum_length
        edge_type = config.edge_type

        if edge_type == EdgeTypes.armchair:
            nanoribbon_length, nanoribbon_width = nanoribbon_width, nanoribbon_length
            vacuum_width, vacuum_length = vacuum_length, vacuum_width

        length_cartesian = nanoribbon_length * np.dot(np.array(material.lattice.vectors[0]), np.array([1, 0, 0]))
        width_cartesian = nanoribbon_width * np.dot(np.array(material.lattice.vectors[1]), np.array([0, 1, 0]))
        height_cartesian = np.dot(np.array(material.lattice.vectors[2]), np.array([0, 0, 1]))
        vacuum_length_cartesian = vacuum_length * np.dot(np.array(material.lattice.vectors[0]), np.array([1, 0, 0]))
        vacuum_width_cartesian = vacuum_width * np.dot(np.array(material.lattice.vectors[1]), np.array([0, 1, 0]))

        return length_cartesian, width_cartesian, height_cartesian, vacuum_length_cartesian, vacuum_width_cartesian

    @staticmethod
    def _get_new_lattice_vectors(
        length: float, width: float, height: float, vacuum_length: float, vacuum_width: float, edge_type: EdgeTypes
    ) -> Tuple[List[float], List[float], List[float]]:
        """
        Calculate the new lattice vectors for the nanoribbon.

        Args:
            length: Length of the nanoribbon.
            width: Width of the nanoribbon.
            height: Height of the nanoribbon.
            vacuum_length: Length of the vacuum region.
            vacuum_width: Width of the vacuum region.
            edge_type: Type of the edge of the nanoribbon.

        Returns:
            Tuple of the new lattice vectors.
        """
        length_lattice_vector = [length + vacuum_length, 0, 0]
        width_lattice_vector = [0, width + vacuum_width, 0]
        height_lattice_vector = [0, 0, height]

        if edge_type == EdgeTypes.armchair:
            length_lattice_vector, width_lattice_vector = width_lattice_vector, length_lattice_vector

        return length_lattice_vector, width_lattice_vector, height_lattice_vector

    @staticmethod
    def _calculate_coordinates_of_cut(
        length: float, width: float, height: float, edge_type: EdgeTypes
    ) -> Tuple[List[float], List[float]]:
        """
        Calculate the coordinates of the rectangular nanoribbon cut from the supercell.

        Args:
            length: Length of the nanoribbon.
            width: Width of the nanoribbon.
            height: Height of the nanoribbon.
            edge_type: Type of the edge of the nanoribbon.

        Returns:
            Tuple of the minimum and maximum coordinates of the cut.
        """
        edge_nudge_value = 0.01
        conditional_nudge_value = edge_nudge_value * (
            -1 * (edge_type == EdgeTypes.armchair) + 1 * (edge_type == EdgeTypes.zigzag)
        )
        min_coordinate = [-edge_nudge_value, conditional_nudge_value, 0]
        max_coordinate = [length - edge_nudge_value, width + conditional_nudge_value, height]
        return min_coordinate, max_coordinate

    def _generate(self, configuration: NanoribbonConfiguration) -> List[_GeneratedItemType]:
        nanoribbon = self.create_nanoribbon(configuration)
        if configuration.edge_type == EdgeTypes.armchair:
            nanoribbon = rotate(nanoribbon, [0, 0, 1], 90)
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
