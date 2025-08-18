from typing import List

from ....utils import get_center_of_coordinates
from ...build_components.entities.reusable.three_dimensional.supercell.helpers import create_supercell
from ...modify import filter_by_condition_on_coordinates
from ...utils import transform_coordinate_to_supercell
from .. import BaseMaterialAnalyzer
from ..coordination import get_voronoi_nearest_neighbors_atom_indices
from ..other import get_closest_site_id_from_coordinate


class CrystalSiteAnalyzer(BaseMaterialAnalyzer):
    coordinate: List[float] = [0.0, 0.0, 0.0]

    @property
    def exact_coordinate(self) -> List[float]:
        return self.coordinate

    @property
    def closest_site_coordinate(self) -> List[float]:
        site_id = get_closest_site_id_from_coordinate(self.material, self.coordinate)
        return self.material.coordinates_array[site_id]

    def get_equidistant_coordinate(self, coordinate=None) -> List[float]:
        """
        Compute a coordinate that is equidistant from the nearest atoms to the target coordinate. Useful for adatom.

        This method works by:
        1. Creating a 3x3x1 supercell of the original material to include atoms in PBC.
        2. Transforming the target coordinate into the supercell's coordinate system.
        3. Using Voronoi tessellation to find the indices of the nearest atoms in the supercell.
        4. Taking the geometric center (average) of these neighboring atoms' coordinates as the equidistant point.
        5. Setting the z-coordinate to the original value (useful for 2D/slab systems).
        6. Transforming the equidistant coordinate back to the original cell's coordinate system.

        The [3, 3, 1] supercell ensures robust neighbor search in x and y, but not in z (for 2D/slab systems).
        """
        if coordinate is None:
            coordinate = self.coordinate
        scaling_factor = [3, 3, 1]
        translation_vector = [1 / 3, 1 / 3, 0]
        supercell_material = create_supercell(self.material, scaling_factor=scaling_factor)

        coordinate_in_supercell = transform_coordinate_to_supercell(
            coordinate=coordinate, scaling_factor=scaling_factor, translation_vector=translation_vector
        )

        neighboring_atoms_ids_in_supercell = get_voronoi_nearest_neighbors_atom_indices(
            material=supercell_material, coordinate=coordinate_in_supercell
        )

        if neighboring_atoms_ids_in_supercell is None:
            raise ValueError("No neighboring atoms found for equidistant calculation.")

        # Filter out atoms that are too close to the z boundaries of the supercell
        supercell_material = filter_by_condition_on_coordinates(supercell_material, lambda c: 1e-2 < c[2] < 1 - 1e-2)
        isolated_neighboring_atoms_basis = supercell_material.basis.model_copy()
        isolated_neighboring_atoms_basis.coordinates.filter_by_ids(neighboring_atoms_ids_in_supercell)
        equidistant_coordinate_in_supercell = get_center_of_coordinates(
            isolated_neighboring_atoms_basis.coordinates.values
        )

        return transform_coordinate_to_supercell(
            equidistant_coordinate_in_supercell, scaling_factor, translation_vector, reverse=True
        )
