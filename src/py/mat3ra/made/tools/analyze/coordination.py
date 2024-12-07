from typing import Dict, List, Optional

import numpy as np
from mat3ra.made.material import Material
from scipy.spatial._ckdtree import cKDTree

from ..convert import to_pymatgen
from ..third_party import PymatgenVoronoiNN
from ..utils import ArrayWithIds, decorator_handle_periodic_boundary_conditions


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


@decorator_handle_periodic_boundary_conditions(0.25)
def get_nearest_neighbors_vectors(
    material: Material,
    indices: Optional[List[int]] = None,
    cutoff: float = 3.0,
    nearest_only: bool = True,
) -> ArrayWithIds:
    """
    Calculate the vectors to the nearest neighbors for each atom in the material.

    Args:
        material (Material): Material object to calculate coordination numbers for.
        indices (List[int]): List of atom indices to calculate coordination numbers for.
        cutoff (float): The maximum cutoff radius for identifying neighbors.
        nearest_only (bool): If True, only consider the first shell of neighbors.

    Returns:
        ArrayWithIds: Array of vectors to the nearest neighbors for each atom.
    """
    new_material = material.clone()
    new_material.to_cartesian()
    if indices is not None:
        new_material.basis.coordinates.filter_by_indices(indices)
    coordinates = np.array(new_material.basis.coordinates.values)
    kd_tree = cKDTree(coordinates)

    nearest_neighbors_vectors = ArrayWithIds()
    for idx, (x, y, z) in enumerate(coordinates):
        # Get all neighbors within cutoff
        distances, neighbors = kd_tree.query(
            [x, y, z], k=len(coordinates), distance_upper_bound=cutoff  # Get all possible neighbors
        )

        if nearest_only:
            # Remove the first distance (distance to self = 0)
            distances = distances[1:]
            neighbors = neighbors[1:]

            # Remove infinite distances (no more neighbors found)
            valid_indices = distances != np.inf
            distances = distances[valid_indices]
            neighbors = neighbors[valid_indices]

            if len(distances) > 0:
                # Find the first shell by analyzing distance distribution
                # Get the first significant gap in distances
                sorted_distances = np.sort(distances)
                distance_gaps = sorted_distances[1:] - sorted_distances[:-1]

                # Find the first significant gap (you might need to adjust this threshold)
                gap_threshold = 0.5  # Å
                significant_gaps = np.where(distance_gaps > gap_threshold)[0]

                if len(significant_gaps) > 0:
                    # Use only neighbors before the first significant gap
                    first_gap_idx = significant_gaps[0] + 1
                    neighbors = neighbors[:first_gap_idx]
        else:
            # Remove self and infinite distances for regular coordination counting
            valid_indices = (distances != np.inf) & (distances != 0)
            neighbors = neighbors[valid_indices]
        vectors = [coordinates[neighbor] - [x, y, z] for neighbor in neighbors]
        nearest_neighbors_vectors.add_item(vectors, idx)

    return nearest_neighbors_vectors


def get_coordination_numbers(material: Material, cutoff: float) -> Dict[int, int]:
    """
    Calculate the coordination numbers for all atoms in the material.

    Args:
        material (Material): The material object.
        cutoff (float): The cutoff radius for identifying neighbors.

    Returns:
        Dict[int, int]: A dictionary mapping atom indices to their coordination numbers.
    """
    nearest_neighbors = get_nearest_neighbors_vectors(material=material, cutoff=cutoff)
    coordination_numbers = {idx: len(vectors) for idx, vectors in enumerate(nearest_neighbors.values)}
    return coordination_numbers


def get_undercoordinated_atom_indices(material: Material, cutoff: float, coordination_threshold: int) -> List[int]:
    """
    Identify undercoordinated atoms based on the coordination threshold (inclusive).

    Args:
        material (Material): The material object.
        cutoff (float): The cutoff radius for identifying neighbors.
        coordination_threshold (int): The coordination number threshold for an atom to be considered undercoordinated.

    Returns:
        List[int]: List of indices of undercoordinated atoms.
    """
    coordination_numbers = get_coordination_numbers(material, cutoff)
    return [idx for idx, number in coordination_numbers.items() if number <= coordination_threshold]


def get_unique_coordination_numbers(material: Material, cutoff: float) -> List[int]:
    """
    Get the unique coordination numbers for all atoms in the material.

    Args:
        material (Material): The material object.
        cutoff (float): The cutoff radius for identifying neighbors.

    Returns:
        List[int]: A list of unique coordination numbers present in the material.
    """
    coordination_numbers = get_coordination_numbers(material, cutoff)
    return sorted(list(set(coordination_numbers.values())))


def are_bonds_templates_similar(template1: np.ndarray, template2: np.ndarray, tolerance: float = 0.1) -> bool:
    """
    Check if two bond templates are similar.

    Args:
        template1 (np.ndarray): First template of bond vectors.
        template2 (np.ndarray): Second template of bond vectors.
        tolerance (float): Angle tolerance for comparison.

    Returns:
        bool: True if the templates are similar, False otherwise.
    """
    if len(template1) != len(template2):
        return False

    dot_matrix = np.dot(template1, template2.T)
    norms1 = np.linalg.norm(template1, axis=1)
    norms2 = np.linalg.norm(template2, axis=1)
    cosine_matrix = dot_matrix / np.outer(norms1, norms2)
    angles_matrix = np.arccos(np.clip(cosine_matrix, -1.0, 1.0))

    unmatched = list(range(len(template2)))
    for angle_row in angles_matrix:
        matches = np.where(angle_row < tolerance)[0]
        if len(matches) == 0:
            return False
        unmatched.remove(matches.tolist()[0])

    return True


def find_template_vectors(material: Material, angle_tolerance: float = 0.1) -> Dict[str, List[np.ndarray]]:
    """
    Find unique bond templates for each element type.

    Args:
        atom_vectors (List[List[np.ndarray]]): List of bond vectors for each atom.
        atom_elements (List[str]): List of chemical elements for each atom.
        angle_tolerance (float): Tolerance for comparing angles between bond vectors.

    Returns:
        Dict[str, List[np.ndarray]]: Dictionary mapping element types to unique bond templates.
    """
    element_templates = {}

    for element in set(material.basis.elements.values):
        element_templates[element] = find_template_vectors_for_element(material, element, angle_tolerance)

    return element_templates


def find_template_vectors_for_element(
    material: Material, element: str, angle_tolerance: float = 0.1
) -> List[np.ndarray]:
    """
    Find unique bond templates for a single element type.

    Args:
        material (Material): The material object.
        element (str): Chemical element to find templates for.
        angle_tolerance (float): Tolerance for comparing angles between bond vectors.

    Returns:
        List[np.ndarray]: List of unique bond templates for the element.
    """
    atom_vectors = get_nearest_neighbors_vectors(material).values
    atom_elements = material.basis.elements.values

    element_indices = [i for i, e in enumerate(atom_elements) if e == element]
    element_vector_lists = [np.array(atom_vectors[i]) for i in element_indices]

    if not element_vector_lists:
        return []

    max_coordination_number = max(len(vectors) for vectors in element_vector_lists)
    max_coordination_number_vectors = [v for v in element_vector_lists if len(v) == max_coordination_number]

    unique_templates: List[np.ndarray] = []
    for template in max_coordination_number_vectors:
        if not any(are_bonds_templates_similar(template, existing, angle_tolerance) for existing in unique_templates):
            unique_templates.append(template)

    return unique_templates


def reconstruct_missing_bonds_for_element(
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


def reconstruct_missing_bonds(
    nearest_neighbor_vectors: List[List[np.ndarray]],
    chemical_elements: List[str],
    templates: Dict[str, List[np.ndarray]],
    angle_tolerance: float = 0.1,
    max_bonds_to_passivate: int = 1,
) -> Dict[int, List[List[float]]]:
    """
    Reconstruct missing bonds for all undercoordinated atoms.

    Args:
        nearest_neighbor_vectors (List[List[np.ndarray]]): List of bond vectors for each atom.
        chemical_elements (List[str]): List of chemical elements for each atom.
        templates (Dict[str, List[np.ndarray]]): Dictionary of bond templates for each element.
        angle_tolerance (float): Tolerance for comparing angles between bond vectors.
        max_bonds_to_passivate (int): Maximum number of bonds to passivate for each undercoordinated atom.

    Returns:
        Dict[int, List[List[float]]]: Dictionary mapping atom indices to reconstructed bond vectors.
    """
    missing_bonds = {}

    for idx, (vectors, element) in enumerate(zip(nearest_neighbor_vectors, chemical_elements)):
        reconstructed_bonds = reconstruct_missing_bonds_for_element(
            vectors, element, templates, angle_tolerance, max_bonds_to_passivate
        )
        if reconstructed_bonds:
            missing_bonds[idx] = reconstructed_bonds

    return missing_bonds
