from typing import List, Optional, Tuple

import numpy as np
from mat3ra.code.array_with_ids import ArrayWithIds
from mat3ra.made.material import Material
from scipy.spatial._ckdtree import cKDTree

from ..bonds import BondDirections, BondDirectionsTemplatesForElement
from ..site import CrystalSite, CrystalSiteList
from .rdf import RadialDistributionFunction
from .utils import decorator_handle_periodic_boundary_conditions


class MaterialWithCrystalSites(Material):
    crystal_sites: CrystalSiteList = CrystalSiteList.from_values([])
    nearest_neighbor_vectors: ArrayWithIds = ArrayWithIds.from_values([])

    @classmethod
    def from_material(cls, material: Material):
        new_material = material.clone()
        new_material.to_cartesian()
        return cls(name=material.name, basis=new_material.basis, lattice=new_material.lattice)

    def analyze(self):
        """
        Analyze the material in place and generate the crystal sites with properties.

        Generated properties:
            - nearest_neighbor_vectors: The nearest neighbor vectors for all sites.
        """
        self.nearest_neighbor_vectors = self.get_neighbors_vectors_for_all_sites(cutoff=3.0)
        self.crystal_sites = CrystalSiteList(
            values=[CrystalSite(nearest_neighbor_vectors=item) for item in self.nearest_neighbor_vectors.values],
            ids=self.nearest_neighbor_vectors.ids,
        )

    @property
    def coordinates_as_kdtree(self):
        return cKDTree(self.basis.coordinates.values)

    def get_neighbors_vectors_for_site(
        self, site_index: int, cutoff: float = 3.0, max_number_of_neighbors: Optional[int] = None
    ) -> List[np.ndarray]:
        """
        Get the nearest neighbor vectors for a specific site.

        Args:
            site_index (int): The index of the site to
            cutoff (float): The maximum cutoff radius for identifying neighbors.
            max_number_of_neighbors (int): The max number of neighbors possible.

        Returns:
            List[np.ndarray]: List of vectors to the nearest neighbors.
        """
        coordinates = self.basis.coordinates.values
        neighbors, distances = self.get_neighbors_for_site(site_index, cutoff, max_number_of_neighbors)
        vectors = [np.array(coordinates[neighbor]) - np.array(coordinates[site_index]) for neighbor in neighbors]
        return vectors

    @decorator_handle_periodic_boundary_conditions(cutoff=0.25)
    def get_neighbors_vectors_for_all_sites(
        self, cutoff: float = 3.0, max_number_of_neighbors: Optional[int] = None
    ) -> ArrayWithIds:
        """
        Get the nearest neighbor vectors for all sites.

        Args:
            cutoff (float): The maximum cutoff radius for identifying neighbors.
            max_number_of_neighbors (int): The max number of neighbors possible.

        Returns:
            ArrayWithIds: Array with the nearest neighbor vectors for all sites.
        """
        nearest_neighbors = ArrayWithIds(values=[], ids=[])
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
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get the nearest neighbors for a specific site.

        Args:
            site_index (int): The index of the site to
            cutoff (float): The maximum cutoff radius for identifying neighbors.
            max_number_of_neighbors (int): The max number of neighbors possible.
            nearest_only (bool): If True, only consider the first shell of neighbors.

        Returns:
            Tuple[np.ndarray, np.ndarray]: Tuple of neighbors and distances.
        """
        coordinates = self.basis.coordinates
        max_number_of_neighbors = (
            len(coordinates.values) if max_number_of_neighbors is None else max_number_of_neighbors
        )
        coordinate = coordinates.get_element_value_by_index(site_index)
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

    def get_nearest_neighbors_for_all_sites(
        self, cutoff: float = 3.0, max_number_of_neighbors: Optional[int] = None
    ) -> ArrayWithIds:
        """
        Get the nearest neighbors for all sites.

        Args:
            cutoff (float): The maximum cutoff radius for identifying neighbors.
            max_number_of_neighbors (int): The max number of neighbors possible.

        Returns:
            ArrayWithIds: Array with the nearest neighbors for all sites.
        """
        nearest_neighbors = ArrayWithIds(values=[], ids=[])
        for site_index in self.basis.coordinates.ids:
            neighbors, distances = self.get_neighbors_for_site(site_index, cutoff, max_number_of_neighbors)
            nearest_neighbors.add_item(neighbors, site_index)
        return nearest_neighbors

    @decorator_handle_periodic_boundary_conditions(cutoff=0.25)
    def get_coordination_numbers(self, cutoff: float = 3.0) -> ArrayWithIds:
        """
        Calculate the coordination numbers for all atoms in the material.

        Args:
            cutoff (float): The cutoff radius for identifying neighbors.

        Returns:
            ArrayWithIds: Array with the coordination numbers for all sites.
        """
        nearest_neighbors = self.get_nearest_neighbors_for_all_sites(cutoff)
        coordination_numbers = [len(neighbors) for neighbors in nearest_neighbors.values]
        return ArrayWithIds(values=coordination_numbers, ids=nearest_neighbors.ids)

    def get_undercoordinated_atom_indices(self, cutoff: float, coordination_threshold: int) -> List[int]:
        """
        Identify undercoordinated atoms based on the coordination threshold (inclusive).

        Args:
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

    @decorator_handle_periodic_boundary_conditions(cutoff=0.25)
    def find_missing_bonds_for_all_sites(
        self,
        bond_templates_list: List[BondDirectionsTemplatesForElement],
        max_bonds_to_add: int = 1,
        angle_tolerance: float = 0.1,
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
        for coordinate in self.basis.coordinates.to_array_of_values_with_ids():
            existing_bond_directions = BondDirections(self.get_neighbors_vectors_for_site(coordinate.id, cutoff=3.0))
            element = self.basis.elements.get_element_value_by_index(coordinate.id)
            bond_templates = [
                bond_templates_for_element.to_ndarray()
                for bond_templates_for_element in bond_templates_list
                if bond_templates_for_element.element == element
            ]
            if bond_templates is None:
                continue
            missing_bonds_for_site = existing_bond_directions.find_missing_directions(
                bond_templates=bond_templates,
                angle_tolerance=angle_tolerance,
                max_bonds_to_add=max_bonds_to_add,
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
            if not any(
                np.array_equal(bond_direction, existing) for existing in unique_templates.bond_directions_templates
            ):
                unique_templates.bond_directions_templates.append(bond_direction)

        return unique_templates
