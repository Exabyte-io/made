import numpy as np

from ......analyze.other import get_surface_atom_indices
from ..enums import SurfaceTypesEnum

from .passivation_material_analyzer import PassivationMaterialAnalyzer


class SurfacePassivationMaterialAnalyzer(PassivationMaterialAnalyzer):
    shadowing_radius: float = 2.5
    depth: float = 5.0

    @property
    def passivant_coordinates(
        self,
    ):
        """
        Calculate the coordinates for placing passivants based on the specified surface type.

        Args:
            self (SurfacePassivationConfiguration): Configuration for passivation.

        Returns:
            list: Coordinates where passivants should be added.
        """
        surface_atoms_indices = get_surface_atom_indices(self.material, self.surface, self.shadowing_radius, self.depth)
        surface_atoms_coordinates = [
            self.material.basis.coordinates.get_element_value_by_index(i) for i in surface_atoms_indices
        ]
        bond_vector = [0, 0, self.bond_length] if self.surface == SurfaceTypesEnum.TOP else [0, 0, -self.bond_length]
        passivant_bond_vector_crystal = self.material.basis.cell.convert_point_to_crystal(bond_vector)
        return (np.array(surface_atoms_coordinates) + np.array(passivant_bond_vector_crystal)).tolist()
