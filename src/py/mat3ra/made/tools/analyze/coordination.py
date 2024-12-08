from enum import Enum
from typing import Dict, List, Optional

import numpy as np
from mat3ra.made.material import Material
from pydantic import BaseModel
from scipy.spatial._ckdtree import cKDTree

from ..convert import to_pymatgen
from ..third_party import PymatgenVoronoiNN
from ..utils import ArrayWithIds, decorator_handle_periodic_boundary_conditions
from . import BaseMaterialAnalyzer


class BondDirectionsTemplatesEnum(List, Enum):
    OCTAHEDRAL = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]]
    TETRAHEDRAL = [[1, 1, 1], [-1, -1, 1], [1, -1, -1], [-1, 1, -1]]
    PLANAR = [[1, 0, 0], [-1, 2, 0], [-1, -2, 0]]
    LINEAR = [[1, 0, 0], [-1, 0, 0]]


class BondDirections(np.ndarray):
    def __eq__(self, other):
        if not isinstance(other, BondDirections):
            return False
        return self.are_bond_directions_similar(self, other)

    @staticmethod
    def are_bond_directions_similar(directions1: np.ndarray, directions2: np.ndarray, tolerance: float = 0.1) -> bool:
        """
        Check if two bond templates are similar.

        Args:
            directions1 (np.ndarray): First template of bond vectors.
            directions2 (np.ndarray): Second template of bond vectors.
            tolerance (float): Angle tolerance for comparison.

        Returns:
            bool: True if the templates are similar, False otherwise.
        """
        if len(directions1) != len(directions2):
            return False

        dot_matrix = np.dot(directions1, directions2.T)
        norms1 = np.linalg.norm(directions1, axis=1)
        norms2 = np.linalg.norm(directions2, axis=1)
        cosine_matrix = dot_matrix / np.outer(norms1, norms2)
        angles_matrix = np.arccos(np.clip(cosine_matrix, -1.0, 1.0))

        unmatched = list(range(len(directions2)))
        for angle_row in angles_matrix:
            matches = np.where(angle_row < tolerance)[0]
            if len(matches) == 0:
                return False
            unmatched.remove(matches.tolist()[0])

        return True

    @staticmethod
    def find_missing_directions(
        vectors: List[np.ndarray],
        element: str,
        templates: Dict[str, List[np.ndarray]],
        angle_tolerance: float = 0.1,
        max_bonds_to_passivate: int = 1,
    ) -> List[List[float]]:
        """
        Reconstruct missing bonds for a single element.

        Args:
            vectors (List[np.ndarray]): List of bond vectors for the atom.
            element (str): Chemical element of the atom.
            templates (Dict[str, List[np.ndarray]]): Dictionary of bond templates for each element.
            angle_tolerance (float): Tolerance for comparing angles between bond vectors.
            max_bonds_to_passivate (int): Maximum number of bonds to passivate for the atom.

        Returns:
            List[List[float]]: List of reconstructed bond vectors.
        """
        if element not in templates:
            return []

        existing_vectors = np.array(vectors) if vectors else np.empty((0, 3))
        max_coordination_number = len(templates[element][0])

        if len(existing_vectors) >= max_coordination_number:
            return []

        best_missing = None
        best_match_count = -1

        for template in templates[element]:
            if existing_vectors.size == 0:
                match_count = 0
            else:
                # TODO: optimize
                dot_matrix = np.dot(template, existing_vectors.T)
                cosine_matrix = dot_matrix / (
                    np.linalg.norm(template, axis=1)[:, None] * np.linalg.norm(existing_vectors, axis=1)
                )
                angles_matrix = np.arccos(np.clip(cosine_matrix, -1.0, 1.0))

                matches = np.any(angles_matrix < angle_tolerance, axis=1)
                match_count = np.sum(matches)

            missing = template[~matches] if existing_vectors.size != 0 else template

            if match_count > best_match_count:
                best_match_count = match_count
                best_missing = missing

        if best_missing is not None:
            num_bonds_to_add = min(
                len(best_missing),
                max_bonds_to_passivate,
                max_coordination_number - len(existing_vectors),
            )
            return best_missing[:num_bonds_to_add].tolist()

        return []


class CrystalSite(BaseModel):
    # element: str
    # coordinate: List[float]
    nearest_neighbor_vectors: List[np.ndarray] = []
    # coordination_number: int = 0
    # see https://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-wp-list for an example
    wyckoff_letter: Optional[str] = None

    class Config:
        arbitrary_types_allowed = True

    @property
    def coordination_number(self):
        return len(self.nearest_neighbor_vectors)


class CrystalSiteList(ArrayWithIds):
    values: List[CrystalSite]


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
            reconstructed_bonds = BondDirections.find_missing_directions(
                vectors, self.basis.elements.values[idx], templates
            )
            if reconstructed_bonds:
                missing_bonds[idx] = reconstructed_bonds
        return missing_bonds


def get_voronoi_nearest_neighbors_atom_indices(
    material: Material,
    coordinate: Optional[List[float]] = None,
    tolerance: float = 0.1,
    cutoff: float = 13.0,
) -> Optional[List[int]]:
    """
    Returns the indices of direct neighboring atoms to a specified position in the material using Voronoi tessellation.

    Args:
        material (Material): The material object to find neighbors in.
        coordinate (List[float]): The position to find neighbors for.
        tolerance (float): tolerance parameter for near-neighbor finding. Faces that are smaller than tol fraction
            of the largest face are not included in the tessellation. (default: 0.1).
            as per: https://pymatgen.org/pymatgen.analysis.html#pymatgen.analysis.local_env.VoronoiNN
        cutoff (float): The cutoff radius for identifying neighbors, in angstroms.

    Returns:
        List[int]: A list of indices of neighboring atoms, or an empty list if no neighbors are found.
    """
    if coordinate is None:
        coordinate = [0, 0, 0]
    structure = to_pymatgen(material)
    voronoi_nn = PymatgenVoronoiNN(
        tol=tolerance,
        cutoff=cutoff,
        weight="solid_angle",
        extra_nn_info=False,
        compute_adj_neighbors=True,
    )
    coordinates = material.basis.coordinates
    site_index = coordinates.get_element_id_by_value(coordinate)
    remove_dummy_atom = False
    if site_index is None:
        structure.append("X", coordinate, validate_proximity=False)
        site_index = len(structure.sites) - 1
        remove_dummy_atom = True
    try:
        neighbors = voronoi_nn.get_nn_info(structure, site_index)
    except ValueError:
        return None
    neighboring_atoms_pymatgen_ids = [n["site_index"] for n in neighbors]
    if remove_dummy_atom:
        structure.remove_sites([-1])

    all_coordinates = material.basis.coordinates
    all_coordinates.filter_by_indices(neighboring_atoms_pymatgen_ids)
    return all_coordinates.ids


def get_undercoordinated_atom_indices(
    material: MaterialWithCrystalSites, cutoff: float, coordination_threshold: int
) -> List[int]:
    """
    Identify undercoordinated atoms based on the coordination threshold (inclusive).

    Args:
        material (MaterialWithCrystalSites): The material object with crystal sites.
        cutoff (float): The cutoff radius for identifying neighbors.
        coordination_threshold (int): The coordination number threshold for an atom to be considered undercoordinated.

    Returns:
        List[int]: List of indices of undercoordinated atoms.
    """
    coordination_numbers = material.get_coordination_numbers(cutoff)
    return [idx for idx, number in coordination_numbers.items() if number <= coordination_threshold]


def get_unique_coordination_numbers(material: MaterialWithCrystalSites, cutoff: float) -> List[int]:
    """
    Get the unique coordination numbers for all atoms in the material.

    Args:
        material (MaterialWithCrystalSites): The material object with crystal sites.
        cutoff (float): The cutoff radius for identifying neighbors.

    Returns:
        List[int]: A list of unique coordination numbers present in the material.
    """
    coordination_numbers = material.get_coordination_numbers(cutoff)
    return sorted(list(set(coordination_numbers.values)))


def find_unique_bond_directions(
    material: MaterialWithCrystalSites, angle_tolerance: float = 0.1
) -> Dict[str, List[np.ndarray]]:
    """
    Find unique bond templates for each element type.

    Args:
        material (Material): The material object.
        angle_tolerance (float): Tolerance for comparing angles between bond vectors.

    Returns:
        Dict[str, List[np.ndarray]]: Dictionary mapping element types to unique bond templates.
    """
    element_templates = {}

    for element in set(material.basis.elements.values):
        element_templates[element] = find_bond_directions_for_element(material, element, angle_tolerance)

    return element_templates


def find_bond_directions_for_element(
    material: MaterialWithCrystalSites, element: str, angle_tolerance: float = 0.1
) -> List[np.ndarray]:
    """
    Find unique bond templates for a single element type.

    Args:
        material (MaterialWithCrystalSites): The material object with crystal sites.
        element (str): Chemical element to find templates for.
        angle_tolerance (float): Tolerance for comparing angles between bond vectors.

    Returns:
        List[np.ndarray]: List of unique bond templates for the element.
    """
    atom_vectors = material.nearest_neighbor_vectors
    atom_elements = material.basis.elements.values

    element_indices = [i for i, e in enumerate(atom_elements) if e == element]
    element_vector_lists = [np.array(atom_vectors[i]) for i in element_indices]

    if not element_vector_lists:
        return []

    max_coordination_number = max(len(vectors) for vectors in element_vector_lists)
    max_coordination_number_vectors = [v for v in element_vector_lists if len(v) == max_coordination_number]

    unique_templates: List[np.ndarray] = []
    for template in max_coordination_number_vectors:
        if all(template != existing for existing in unique_templates):
            unique_templates.append(template)

    return unique_templates


####################################################################################################
# Radial Distribution Function (RDF) Analysis
####################################################################################################


class RadicalDistributionFunction(BaseModel):
    rdf: np.ndarray
    bin_centers: np.ndarray

    class Config:
        arbitrary_types_allowed = True

    @classmethod
    def from_material(cls, material: Material, cutoff: float = 10.0, bin_size: float = 0.1):
        analyzer = BaseMaterialAnalyzer(material)
        distances = analyzer.pairwise_distances
        density = analyzer.atomic_density
        distances = distances[distances <= cutoff]

        # Bin distances into a histogram
        bins = np.arange(0, cutoff + bin_size, bin_size)  # Bin edges
        hist, bin_edges = np.histogram(distances, bins=bins, density=False)

        # Convert to radial distribution function
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        shell_volumes = (
            (4 / 3) * np.pi * (np.power(bin_edges[1:], 3) - np.power(bin_edges[:-1], 3))
        )  # Volume of spherical shells

        rdf = hist / (shell_volumes * density)

        return cls(rdf=rdf, bin_centers=bin_centers)

    @property
    def first_peak_index(self):
        return np.argmax(self.rdf[1:]) + 1

    @property
    def first_peak_value(self):
        return self.rdf[self.first_peak_index]

    @property
    def first_peak_width(self):
        half_max = 0.5 * self.first_peak_value

        # Find left boundary
        left_index = self.first_peak_index
        while left_index > 0 and self.rdf[left_index] > half_max:
            left_index -= 1

        # Find right boundary
        right_index = self.first_peak_index
        while right_index < len(self.rdf) - 1 and self.rdf[right_index] > half_max:
            right_index += 1

        # Compute the width
        left_boundary = self.bin_centers[left_index]
        right_boundary = self.bin_centers[right_index]
        first_peak_width = right_boundary - left_boundary

        return float(first_peak_width)

    @property
    def first_peak_distance(self):
        return self.bin_centers[self.first_peak_index]

    def is_within_first_peak(self, distance: float, tolerance: float = 0.1) -> bool:
        return (
            self.first_peak_distance - 0.5 * self.first_peak_width - tolerance
            < distance
            < self.first_peak_distance + 0.5 * self.first_peak_width + tolerance
        )
