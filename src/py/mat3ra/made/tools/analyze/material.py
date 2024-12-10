from typing import List, Optional

import numpy as np
from mat3ra.made.material import Material
from mat3ra.made.utils import ArrayWithIds
from scipy.spatial._ckdtree import cKDTree

from ..bonds import BondDirections, BondDirectionsTemplatesForElement
from ..site import CrystalSite, CrystalSiteList
from .rdf import RadialDistributionFunction
from .utils import decorator_handle_periodic_boundary_conditions


class MaterialWithCrystalSites(Material):
    crystal_sites: CrystalSiteList = CrystalSiteList(values=[])

    @classmethod
    def from_material(cls, material: Material):
        material.to_cartesian()
        config = material.to_json()
        return cls(config)

    def analyze(self):
        self.nearest_neighbor_vectors = self.get_neighbors_vectors_for_all_sites(cutoff=3.0)
        self.crystal_sites = CrystalSiteList(
            values=[CrystalSite(nearest_neighbor_vectors=item) for item in self.nearest_neighbor_vectors.values],
            ids=self.nearest_neighbor_vectors.ids,
        )

    @property
    def coordinates_as_kdtree(self):
        return cKDTree(self.basis.coordinates.values)

    @decorator_handle_periodic_boundary_conditions(cutoff=0.25)
    def get_neighbors_vectors_for_all_sites(
        self, cutoff: float = 3.0, max_number_of_neighbors: Optional[int] = None
    ) -> ArrayWithIds:
        nearest_neighbors = ArrayWithIds()
        for site_index in range(len(self.basis.coordinates.values)):
            vectors = self.get_neighbors_vectors_for_site(site_index, cutoff, max_number_of_neighbors)
            nearest_neighbors.add_item(vectors, site_index)
        return nearest_neighbors

    def get_neighbors_vectors_for_site(
        self, site_index: int, cutoff: float = 3.0, max_number_of_neighbors: Optional[int] = None
    ) -> List[np.ndarray]:
        coordinates = self.basis.coordinates.values
        neighbors, distances = self.get_neighbors_for_site(site_index, cutoff, max_number_of_neighbors)
        vectors = [np.array(coordinates[neighbor]) - np.array(coordinates[site_index]) for neighbor in neighbors]
        return vectors

    def get_nearest_neighbors_for_all_sites(
        self, cutoff: float = 3.0, max_number_of_neighbors: Optional[int] = None
    ) -> ArrayWithIds:
        """
        Get the nearest neighbors for all sites.

        Args:
            cutoff (float): The maximum cutoff radius for identifying neighbors.
            max_number_of_neighbors (int): The max number of neighbors possible.
        """
        nearest_neighbors = ArrayWithIds()
        for site_index in range(len(self.basis.coordinates.values)):
            neighbors, distances = self.get_neighbors_for_site(site_index, cutoff, max_number_of_neighbors)
            nearest_neighbors.add_item(neighbors, site_index)
        return nearest_neighbors

    def get_neighbors_for_site(
        self,
        site_index: int,
        cutoff: float = 3.0,
        max_number_of_neighbors: Optional[int] = None,
        nearest_only: bool = True,
    ):
        """
        Get the nearest neighbors for a specific site.

        Args:
            site_index (int): The index of the site to
            cutoff (float): The maximum cutoff radius for identifying neighbors.
            max_number_of_neighbors (int): The max number of neighbors possible.
            nearest_only (bool): If True, only consider the first shell of neighbors.
            distance_tolerance (float): The tolerance for identifying nearest_neighbors.
        """
        coordinates = self.basis.coordinates.values
        max_number_of_neighbors = len(coordinates) if max_number_of_neighbors is None else max_number_of_neighbors
        coordinate = coordinates[site_index]
        distances, neighbors = self.coordinates_as_kdtree.query(
            coordinate, k=max_number_of_neighbors, distance_upper_bound=cutoff
        )

        valid_indices = (distances != np.inf) & (distances != 0)
        distances = distances[valid_indices]
        neighbors = neighbors[valid_indices]

        if nearest_only:
            rdf = RadialDistributionFunction.from_material(self)
            valid_indices = np.where([rdf.is_within_first_peak(distance) for distance in distances])
            distances = distances[valid_indices]
            neighbors = neighbors[valid_indices]
        return neighbors, distances

    @decorator_handle_periodic_boundary_conditions(cutoff=0.25)
    def get_coordination_numbers(self, cutoff: float = 3.0) -> ArrayWithIds:
        """
        Calculate the coordination numbers for all atoms in the material.

        Args:
            cutoff (float): The cutoff radius for identifying neighbors.

        Returns:
            Dict[int, int]: A dictionary mapping atom indices to their coordination numbers.
        """
        nearest_neighbors = self.get_nearest_neighbors_for_all_sites(cutoff)
        coordination_numbers = [len(neighbors) for neighbors in nearest_neighbors.values]
        return ArrayWithIds(values=coordination_numbers, ids=nearest_neighbors.ids)

    def get_undercoordinated_atom_indices(self, cutoff: float, coordination_threshold: int) -> List[int]:
        """
        Identify undercoordinated atoms based on the coordination threshold (inclusive).

        Args:
            material (MaterialWithCrystalSites): The material object with crystal sites.
            cutoff (float): The cutoff radius for identifying neighbors.
            coordination_threshold (int): The coordination number threshold for an atom
                to be considered undercoordinated, inclusive.

        Returns:
            List[int]: List of indices of undercoordinated atoms.
        """
        coordination_numbers = self.get_coordination_numbers(cutoff)
        return [
            item.id
            for item in coordination_numbers.to_array_of_values_with_ids()
            if item.value <= coordination_threshold
        ]

    def get_unique_coordination_numbers(self, cutoff: float) -> List[int]:
        """
        Get the unique coordination numbers for all atoms in the material.

        Args:
            self (MaterialWithCrystalSites): The material object with crystal sites.
            cutoff (float): The cutoff radius for identifying neighbors.

        Returns:
            List[int]: A list of unique coordination numbers present in the material.
        """
        coordination_numbers = self.get_coordination_numbers(cutoff)
        return sorted(list(set(coordination_numbers.values)))

    def find_missing_bonds_for_all_sites(
        self, bond_templates_list: List[BondDirectionsTemplatesForElement]
    ) -> List[BondDirections]:
        """
        Find missing bonds for all sites in the material.

        Args:
            self (MaterialWithCrystalSites): The material object with crystal sites.
            bond_templates_list (List[BondDirectionsTemplatesForElement]): List of bond templates for each element.

        Returns:
            List[BondDirections]: List of missing bonds for each site.
        """
        missing_bonds: List[BondDirections] = []
        for id, coordinate in self.basis.coordinates.to_array_of_values_with_ids():
            existing_bond_directions = BondDirections(self.get_neighbors_vectors_for_site(id, cutoff=3.0))
            bond_templates = next((b for b in bond_templates_list if b.element == self.basis.elements.values[id]), None)
            if bond_templates is None:
                continue

            missing_bonds_for_site = existing_bond_directions.find_missing_directions(
                templates=bond_templates.to_ndarray(), angle_tolerance=0.1, max_bonds_to_add=1
            )
            missing_bonds.append(BondDirections(missing_bonds_for_site))

        return missing_bonds

    def find_unique_bond_directions(self) -> List[BondDirectionsTemplatesForElement]:
        """
        Find unique bond templates for each element type.

        Args:
            self (Material): The material object.

        Returns:
            List[BondDirectionsTemplatesForElement]: List of unique bond templates for each element.

        Example schematic structure for Graphene with 2 unique bond directions:
            {"C": [[0, 1, 0], [-2, -1, 0], [2, -1, 0]], [[0,-1,0], [-2, 1, 0], [2, 1, 0]]}
        """
        element_templates: List[BondDirectionsTemplatesForElement] = []

        for element in set(self.basis.elements.values):
            element_templates.append(self.find_bond_directions_for_element(element))

        return element_templates

    def find_bond_directions_for_element(self, element: str) -> BondDirectionsTemplatesForElement:
        """
        Find bond directions templates for each element.

        Args:
            self (MaterialWithCrystalSites): The material object with crystal sites.
            element (str): Chemical element to find templates for.

        Returns:
            List[BondDirectionsTemplatesForElement]: List of unique bond templates for the element.
        """
        atom_vectors = self.nearest_neighbor_vectors
        atom_elements = self.basis.elements.values

        element_indices = [i for i, e in enumerate(atom_elements) if e == element]
        element_vector_lists = [np.array(atom_vectors.get_element_value_by_index(i)) for i in element_indices]

        if not element_vector_lists:
            return BondDirectionsTemplatesForElement(bond_directions_templates=[], element=element)

        max_coordination_number = max(len(vectors) for vectors in element_vector_lists)
        max_coordination_number_vectors = [v for v in element_vector_lists if len(v) == max_coordination_number]

        unique_templates: BondDirectionsTemplatesForElement = BondDirectionsTemplatesForElement(
            bond_directions_templates=[], element=element
        )
        for template in max_coordination_number_vectors:
            bond_direction = BondDirections(template)
            if all(bond_direction != existing for existing in unique_templates):
                unique_templates.bond_directions_templates.append(bond_direction)

        return unique_templates
