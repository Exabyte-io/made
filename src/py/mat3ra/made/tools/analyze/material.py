from typing import Dict, List, Optional

import numpy as np
from mat3ra.made.material import Material
from mat3ra.made.utils import ArrayWithIds
from scipy.spatial._ckdtree import cKDTree

from ..bonds import BondDirections
from ..site import CrystalSite, CrystalSiteList
from .rdf import RadicalDistributionFunction
from .utils import decorator_handle_periodic_boundary_conditions


class MaterialWithCrystalSites(Material):
    crystal_sites: CrystalSiteList = CrystalSiteList(values=[])

    def __init__(self, **data):
        super().__init__(config=data)
        self.nearest_neighbor_vectors = self.get_neighbors_vectors_for_all_sites(cutoff=3.0)
        self.crystal_sites = CrystalSiteList(
            values=[CrystalSite(nearest_neighbor_vectors=item) for item in self.nearest_neighbor_vectors.values],
            ids=self.nearest_neighbor_vectors.ids,
        )

    @classmethod
    def from_material(cls, material: Material):
        config = material.to_json()
        return cls(**config)

    @property
    def coordinates_as_kdtree(self):
        return cKDTree(self.basis.coordinates.values)

    # @decorator_handle_periodic_boundary_conditions(cutoff=0.25)
    def get_neighbors_vectors_for_site(
        self, site_index: int, cutoff: float = 3.0, max_number_of_neighbors: Optional[int] = None
    ) -> List[np.ndarray]:
        coordinates = self.basis.coordinates.values
        neighbors, distances = self.get_neighbors_for_site(site_index, cutoff, max_number_of_neighbors)
        vectors = [np.array(coordinates[neighbor]) - np.array(coordinates[site_index]) for neighbor in neighbors]
        return vectors

    def get_neighbors_vectors_for_all_sites(
        self, cutoff: float = 3.0, max_number_of_neighbors: Optional[int] = None
    ) -> ArrayWithIds:
        nearest_neighbors = ArrayWithIds()
        for site_index in range(len(self.basis.coordinates.values)):
            vectors = self.get_neighbors_vectors_for_site(site_index, cutoff, max_number_of_neighbors)
            nearest_neighbors.add_item(vectors, site_index)
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
            rdf = RadicalDistributionFunction.from_material(self)
            valid_indices = np.where([rdf.is_within_first_peak(distance) for distance in distances])
            distances = distances[valid_indices]
            neighbors = neighbors[valid_indices]
        return neighbors, distances

    @decorator_handle_periodic_boundary_conditions(0.25)
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

    def get_coordination_numbers(self, cutoff: float = 3.0) -> ArrayWithIds:
        """
        Calculate the coordination numbers for all atoms in the material.

        Args:
            cutoff (float): The cutoff radius for identifying neighbors.

        Returns:
            Dict[int, int]: A dictionary mapping atom indices to their coordination numbers.
        """
        nearest_neighbors = self.get_nearest_neighbors_for_all_sites(cutoff)
        coordination_numbers = ArrayWithIds(values=nearest_neighbors)
        return coordination_numbers

    def find_missing_bonds_for_all_sites(self, templates: Dict[str, List[np.ndarray]]) -> Dict[int, List[List[float]]]:
        missing_bonds = {}
        for idx, (vectors) in enumerate(self.nearest_neighbor_vectors):
            reconstructed_bonds = BondDirections.find_missing_directions(vectors, templates)
            if reconstructed_bonds:
                missing_bonds[idx] = reconstructed_bonds
        return missing_bonds

    def get_undercoordinated_atom_indices(self, cutoff: float, coordination_threshold: int) -> List[int]:
        """
        Identify undercoordinated atoms based on the coordination threshold (inclusive).

        Args:
            material (MaterialWithCrystalSites): The material object with crystal sites.
            cutoff (float): The cutoff radius for identifying neighbors.
            coordination_threshold (int): The coordination number threshold for an atom to be considered undercoordinated.

        Returns:
            List[int]: List of indices of undercoordinated atoms.
        """
        coordination_numbers = self.get_coordination_numbers(cutoff)
        return [idx for idx, number in coordination_numbers.items() if number <= coordination_threshold]

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

    def find_unique_bond_directions(self, angle_tolerance: float = 0.1) -> Dict[str, List[BondDirections]]:
        """
        Find unique bond templates for each element type.

        Args:
            self (Material): The material object.
            angle_tolerance (float): Tolerance for comparing angles between bond vectors.

        Returns:
            Dict[str, List[np.ndarray]]: Dictionary mapping element types to unique bond templates.
        """
        element_templates = {}

        for element in set(self.basis.elements.values):
            element_templates[element] = self.find_bond_directions_for_element(element, angle_tolerance)

        return element_templates

    def find_bond_directions_for_element(self, element: str, angle_tolerance: float = 0.1) -> List[BondDirections]:
        """
        Find unique bond templates for a single element type.

        Args:
            self (MaterialWithCrystalSites): The material object with crystal sites.
            element (str): Chemical element to find templates for.
            angle_tolerance (float): Tolerance for comparing angles between bond vectors.

        Returns:
            List[BondDirections]: List of unique bond templates for the element.
        """
        atom_vectors = self.nearest_neighbor_vectors
        atom_elements = self.basis.elements.values

        element_indices = [i for i, e in enumerate(atom_elements) if e == element]
        element_vector_lists = [np.array(atom_vectors[i]) for i in element_indices]

        if not element_vector_lists:
            return []

        max_coordination_number = max(len(vectors) for vectors in element_vector_lists)
        max_coordination_number_vectors = [v for v in element_vector_lists if len(v) == max_coordination_number]

        unique_templates: List[BondDirections] = []
        for template in max_coordination_number_vectors:
            bond_direction = BondDirections(template)
            if all(bond_direction != existing for existing in unique_templates):
                unique_templates.append(bond_direction)

        return unique_templates
