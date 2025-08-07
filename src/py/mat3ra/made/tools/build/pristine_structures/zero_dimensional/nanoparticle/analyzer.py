from typing import List, Union

import numpy as np

from .....analyze import BaseMaterialAnalyzer
from .....entities.coordinate import SphereCoordinateCondition, CoordinateCondition


class NanoparticleMaterialAnalyzer(BaseMaterialAnalyzer):
    orientation_z: List[int] = [0, 0, 1]
    vacuum_padding: float = 10.0
    coordinate_condition: Union[SphereCoordinateCondition, CoordinateCondition] = SphereCoordinateCondition()

    def _calculate_max_dimension_from_condition(self) -> float:
        if isinstance(self.coordinate_condition, SphereCoordinateCondition):
            return 2.0 * self.coordinate_condition.radius

        return self.vacuum_padding

    def _get_cartesian_lattice_dimensions(self) -> List[float]:
        """
        Get the Cartesian dimensions of the unit cell lattice vectors.
        """
        lattice_vectors = self.material.lattice.vector_arrays
        return [np.linalg.norm(vector) for vector in lattice_vectors]

    @property
    def supercell_to_fit_nanoparticle(self) -> List[List[int]]:
        """
        Calculate the supercell size needed to fit the nanoparticle.
        This calculates how many unit cells are needed in each direction to contain
        the nanoparticle defined by the coordinate condition.
        """
        max_dimension = self._calculate_max_dimension_from_condition()

        lattice_dimensions = self._get_cartesian_lattice_dimensions()

        nx = max(1, int(np.ceil(max_dimension / lattice_dimensions[0])))
        ny = max(1, int(np.ceil(max_dimension / lattice_dimensions[1])))
        nz = max(1, int(np.ceil(max_dimension / lattice_dimensions[2])))

        return [[nx, 0, 0], [0, ny, 0], [0, 0, nz]]

    @property
    def supercell_for_vacuum(self) -> List[List[int]]:
        """
        Calculate the supercell size that includes vacuum padding on all sides.
        This is used to create the final cell with appropriate vacuum space.
        """
        nanoparticle_supercell = self.supercell_to_fit_nanoparticle

        # Get base dimensions from nanoparticle supercell
        nx_base = nanoparticle_supercell[0][0]
        ny_base = nanoparticle_supercell[1][1]
        nz_base = nanoparticle_supercell[2][2]

        # Get lattice dimensions
        lattice_dimensions = self._get_cartesian_lattice_dimensions()

        # Calculate additional cells needed for vacuum padding on each side
        # We need padding on both sides, so total padding = 2 * vacuum_padding
        vacuum_cells_x = max(1, int(np.ceil(2 * self.vacuum_padding / lattice_dimensions[0])))
        vacuum_cells_y = max(1, int(np.ceil(2 * self.vacuum_padding / lattice_dimensions[1])))
        vacuum_cells_z = max(1, int(np.ceil(2 * self.vacuum_padding / lattice_dimensions[2])))

        # Total supercell size = nanoparticle size + vacuum padding
        nx_total = nx_base + vacuum_cells_x
        ny_total = ny_base + vacuum_cells_y
        nz_total = nz_base + vacuum_cells_z

        return [[nx_total, 0, 0], [0, ny_total, 0], [0, 0, nz_total]]

    @property
    def xy_supercell_matrix_to_fit_nanoparticle(self) -> List[List[int]]:
        """
        Get the 2D supercell matrix for the xy plane to fit the nanoparticle.
        """
        supercell_3d = self.supercell_to_fit_nanoparticle
        return [[supercell_3d[0][0], 0], [0, supercell_3d[1][1]]]

    @property
    def number_of_layers_to_fit_nanoparticle(self) -> int:
        """
        Get the number of layers (z-direction repetitions) needed to fit the nanoparticle.
        """
        return self.supercell_to_fit_nanoparticle[2][2]
