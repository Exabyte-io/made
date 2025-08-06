from typing import List

import numpy as np

from ......analyze.material import MaterialWithCrystalSites
from ......bond_directions.bond_directions import BondDirections
from ......bond_directions.bond_directions_templates_for_element import BondDirectionsTemplatesForElement

from .surface_passivation_material_analyzer import SurfacePassivationMaterialAnalyzer


class CoordinationBasedPassivationMaterialAnalyzer(SurfacePassivationMaterialAnalyzer):
    """
    Builder for passivating material surfaces based on coordination number analysis.

        This builder analyzes atomic coordination environments and reconstructs missing bonds
        using templates derived from fully coordinated atoms. It works by:

        1. Finding undercoordinated atoms that need passivation
        2. Creating bond vector templates for each element type by:
            - Collecting vectors from atoms with the highest valid coordination
            - Grouping similar vector arrangements into unique templates
        3. Reconstructing missing bonds by:
            - Matching existing bonds against templates
            - Finding the best-matching template for each atom
            - Adding missing bond vectors from the template
        4. Placing passivant atoms (typically H) along the reconstructed bond vectors
    """

    coordination_threshold: int = 3  # The coordination number threshold for an atom to be considered undercoordinated.
    number_of_bonds_to_passivate: int = 1  # The maximum number of bonds to passivate for each undercoordinated atom.
    symmetry_tolerance: float = 0.1  # The tolerance for symmetry comparison of vectors for bonds.

    @property
    def material_with_crystal_sites(self) -> MaterialWithCrystalSites:
        material_with_crystal_sites = MaterialWithCrystalSites.from_material(material=self.material)
        material_with_crystal_sites.analyze()
        return material_with_crystal_sites

    @property
    def undercoordinated_atoms_indices(self):
        return self.material_with_crystal_sites.get_undercoordinated_atom_indices(
            cutoff=self.shadowing_radius,
            coordination_threshold=self.coordination_threshold,
        )

    @property
    def bonds_templates(self) -> List[BondDirectionsTemplatesForElement]:
        """
        Find unique bond directions based on the material's crystal sites.

        Returns:
            List[BondDirectionsTemplatesForElement]: A list of unique bond direction templates.
        """
        # TODO: bonds_templates will be passed from the self in the "controlled" version of this class
        return self.material_with_crystal_sites.find_unique_bond_directions()

    @property
    def reconstructed_bonds(self) -> List[BondDirections]:
        """
        Find missing bonds for all sites based on the bond templates.

        Returns:
            List[BondDirections]: A list of reconstructed bond directions for each atom.
        """
        return self.material_with_crystal_sites.find_missing_bonds_for_all_sites(
            bond_templates_list=self.bonds_templates,
            max_bonds_to_add=self.number_of_bonds_to_passivate,
            angle_tolerance=self.symmetry_tolerance,
        )

    def _normalize_and_convert_bond_vector(self, bond_vector):
        """
        Normalize a bond vector and convert it to crystal coordinates.
        Args:
            bond_vector: The bond vector to normalize and convert.

        Returns:
            numpy.ndarray: The normalized bond vector in crystal coordinates, or None if the vector is zero.
        """
        bond_vector_np = np.array(bond_vector)
        if np.linalg.norm(bond_vector_np) == 0:
            return None
        normalized_bond = bond_vector_np / np.linalg.norm(bond_vector_np) * self.bond_length
        return self.material.basis.cell.convert_point_to_crystal(normalized_bond)

    @property
    def passivant_coordinates(
        self,
    ):
        """
        Calculate the coordinates for placing passivating atoms based on reconstructed bonds.

        Args:
            self (PassivationConfiguration): Configuration for passivation.

        Returns:
            list: Coordinates where passivants should be added.
        """
        passivant_coordinates = []

        for idx in self.undercoordinated_atoms_indices:
            for bond_vector in self.reconstructed_bonds[idx]:
                normalized_bond_crystal = self._normalize_and_convert_bond_vector(bond_vector)
                if normalized_bond_crystal is None:
                    continue
                passivant_coordinates.append(
                    np.array(self.material.basis.coordinates.get_element_value_by_index(idx))
                    + np.array(normalized_bond_crystal)
                )

        return passivant_coordinates
