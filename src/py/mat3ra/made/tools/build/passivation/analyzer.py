from typing import List

import numpy as np

from ...analyze.material import MaterialWithCrystalSites
from ...analyze.other import get_surface_atom_indices
from ...analyze.slab import SlabMaterialAnalyzer
from ...bonds import BondDirections
from ...enums import SurfaceTypes


class PassivationMaterialAnalyzer(SlabMaterialAnalyzer):
    passivant: str = "H"
    bond_length: float = 1.0
    surface: SurfaceTypes = SurfaceTypes.TOP

    def get_passivant_coordinates(self):
        raise NotImplementedError(
            "This method should be implemented in subclasses to return passivant coordinates based on the surface type."
        )


class SurfacePassivationMaterialAnalyzer(PassivationMaterialAnalyzer):
    shadowing_radius: float = 2.5
    depth: float = 5.0

    def get_passivant_coordinates(
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
        bond_vector = [0, 0, self.bond_length] if self.surface == SurfaceTypes.TOP else [0, 0, -self.bond_length]
        passivant_bond_vector_crystal = self.material.basis.cell.convert_point_to_crystal(bond_vector)
        return (np.array(surface_atoms_coordinates) + np.array(passivant_bond_vector_crystal)).tolist()


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
        return MaterialWithCrystalSites.from_material(material=self.material)

    def get_passivant_coordinates(
        self,
    ):
        """
        Calculate the coordinates for placing passivating atoms based on reconstructed bonds.

        Args:
            self (PassivationConfiguration): Configuration for passivation.

        Returns:
            list: Coordinates where passivants should be added.
        """
        material_with_crystal_sites = self.material_with_crystal_sites
        material_with_crystal_sites.analyze()

        undercoordinated_atoms_indices = material_with_crystal_sites.get_undercoordinated_atom_indices(
            cutoff=self.shadowing_radius,
            coordination_threshold=self.coordination_threshold,
        )
        # TODO: bonds_templates will be passed from the self in the "controlled" version of this class
        bonds_templates = material_with_crystal_sites.find_unique_bond_directions()
        reconstructed_bonds: List[BondDirections] = material_with_crystal_sites.find_missing_bonds_for_all_sites(
            bond_templates_list=bonds_templates,
            max_bonds_to_add=self.number_of_bonds_to_passivate,
            angle_tolerance=self.symmetry_tolerance,
        )

        passivant_coordinates = []

        for idx in undercoordinated_atoms_indices:
            for bond_vector in reconstructed_bonds[idx]:
                bond_vector_np = np.array(bond_vector)
                if np.linalg.norm(bond_vector_np) == 0:
                    continue  # Avoid division by zero
                normalized_bond = bond_vector_np / np.linalg.norm(bond_vector_np) * self.bond_length
                normalized_bond_crystal = self.material.basis.cell.convert_point_to_crystal(normalized_bond)
                passivant_coordinates.append(
                    np.array(self.material.basis.coordinates.get_element_value_by_index(idx))
                    + np.array(normalized_bond_crystal)
                )

        return passivant_coordinates
