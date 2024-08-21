from typing import Callable, List, Literal, Optional, Tuple

import numpy as np

from ..material import Material
from .convert import decorator_convert_material_args_kwargs_to_atoms, to_pymatgen
from .third_party import ASEAtoms, PymatgenIStructure, PymatgenVoronoiNN


@decorator_convert_material_args_kwargs_to_atoms
def get_average_interlayer_distance(
    interface_atoms: ASEAtoms, tag_substrate: str, tag_film: str, threshold: float = 0.5
) -> float:
    """
    Calculate the average distance between the top layer of substrate atoms and the bottom layer of film atoms.

    Args:
        interface_atoms (ase.ASEAtoms): The ASE ASEAtoms object containing both sets of atoms.
        tag_substrate (int): The tag representing the substrate atoms.
        tag_film (int): The tag representing the film atoms.
        threshold (float): The threshold for identifying the top and bottom layers of atoms.

    Returns:
        float: The average distance between the top layer of substrate and the bottom layer of film.
    """
    # Extract z-coordinates of substrate and film atoms
    z_substrate = interface_atoms.positions[interface_atoms.get_tags() == tag_substrate][:, 2]
    z_film = interface_atoms.positions[interface_atoms.get_tags() == tag_film][:, 2]

    # Identify the top layer of substrate atoms and bottom layer of film atoms by z-coordinate
    top_substrate_layer_z = np.max(z_substrate)
    bottom_film_layer_z = np.min(z_film)

    # Get the average z-coordinate of the top substrate atoms (within a threshold from the top)
    top_substrate_atoms = z_substrate[z_substrate >= top_substrate_layer_z - threshold]
    avg_z_top_substrate = np.mean(top_substrate_atoms)

    # Get the average z-coordinate of the bottom film atoms (within a threshold from the bottom)
    bottom_film_atoms = z_film[z_film <= bottom_film_layer_z + threshold]
    avg_z_bottom_film = np.mean(bottom_film_atoms)

    # Calculate the average distance between the top layer of substrate and the bottom layer of film
    average_interlayer_distance = avg_z_bottom_film - avg_z_top_substrate
    return abs(average_interlayer_distance)


@decorator_convert_material_args_kwargs_to_atoms
def get_surface_area(atoms: ASEAtoms):
    """
    Calculate the area of the surface perpendicular to the z-axis of the atoms structure.

    Args:
        atoms (ase.ASEAtoms): The ASEAtoms object to calculate the surface area of.

    Returns:
        float: The surface area of the atoms.
    """
    matrix = atoms.cell
    cross_product = np.cross(matrix[0], matrix[1])
    return np.linalg.norm(cross_product)


@decorator_convert_material_args_kwargs_to_atoms
def get_chemical_formula(atoms: ASEAtoms):
    """
    Calculate the formula of the atoms structure.

    Args:
        atoms (ase.ASEAtoms): The ASEAtoms object to calculate the formula of.

    Returns:
        str: The formula of the atoms.
    """
    return atoms.get_chemical_formula()


def get_closest_site_id_from_coordinate(material: Material, coordinate: List[float]) -> int:
    """
    Get the site ID of the closest site to a given coordinate in the crystal.

    Args:
        material (Material): The material object to find the closest site in.
        coordinate (List[float]): The coordinate to find the closest site to.

    Returns:
        int: The site ID of the closest site.
    """
    coordinates = np.array(material.coordinates_array)
    coordinate = np.array(coordinate)  # type: ignore
    distances = np.linalg.norm(coordinates - coordinate, axis=1)
    return int(np.argmin(distances))


def get_closest_site_id_from_coordinate_and_element(
    material: Material, coordinate: List[float], chemical_element: Optional[str] = None
) -> int:
    """
    Get the site ID of the closest site with a given element to a given coordinate in the crystal.

    Args:
        material (Material): The material object to find the closest site in.
        coordinate (List[float]): The coordinate to find the closest site to.
        chemical_element (str): The element of the site to find.

    Returns:
        int: The site ID of the closest site with the given element.
    """
    coordinates = np.array(material.basis.coordinates.values)
    coordinate = np.array(coordinate)  # type: ignore
    distances = np.linalg.norm(coordinates - coordinate, axis=1)

    if chemical_element is not None:
        elements = np.array(material.basis.elements.values)
        element_indices = np.where(elements == chemical_element)[0]
        distances = distances[element_indices]
        return int(element_indices[np.argmin(distances)])
    else:
        return int(np.argmin(distances))


def get_atom_indices_within_layer_by_atom_index(material: Material, atom_index: int, layer_thickness: float):
    """
    Select all atoms within a specified layer thickness of a central atom along a direction.
    This direction will be orthogonal to the AB plane.
    Layer thickness is converted from angstroms to fractional units based on the lattice vector length.

    Args:
        material (Material): Material object
        atom_index (int): Index of the central atom
        layer_thickness (float): Thickness of the layer in angstroms

    Returns:
        List[int]: List of indices of atoms within the specified layer
    """
    coordinates = material.basis.coordinates.to_array_of_values_with_ids()
    vectors = material.lattice.vectors
    direction_vector = np.array(vectors[2])

    # Normalize the direction vector
    direction_length = np.linalg.norm(direction_vector)
    direction_norm = direction_vector / direction_length
    central_atom_coordinate = coordinates[atom_index]
    central_atom_projection = np.dot(central_atom_coordinate.value, direction_norm)

    layer_thickness_frac = layer_thickness / direction_length

    lower_bound = central_atom_projection - layer_thickness_frac / 2
    upper_bound = central_atom_projection + layer_thickness_frac / 2

    selected_indices = []
    for coord in coordinates:
        # Project each coordinate onto the direction vector
        projection = np.dot(coord.value, direction_norm)
        if lower_bound <= projection <= upper_bound:
            selected_indices.append(coord.id)
    return selected_indices


def get_atom_indices_within_layer_by_atom_coordinate(
    material: Material, coordinate: List[float], layer_thickness: float
):
    """
    Select all atoms within a specified layer thickness of a central atom along a direction.
    This direction will be orthogonal to the AB plane.
    Layer thickness is converted from angstroms to fractional units based on the lattice vector length.

    Args:
        material (Material): Material object
        coordinate (List[float]): Coordinate of the central atom in crystal coordinates
        layer_thickness (float): Thickness of the layer in angstroms

    Returns:
        List[int]: List of indices of atoms within the specified layer
    """
    site_id = get_closest_site_id_from_coordinate(material, coordinate)
    return get_atom_indices_within_layer_by_atom_index(material, site_id, layer_thickness)


def get_atom_indices_within_layer(
    material: Material,
    atom_index: Optional[int] = 0,
    coordinate: Optional[List[float]] = None,
    layer_thickness: float = 1,
):
    """
    Select all atoms within a specified layer thickness of the central atom along the c-vector direction.

    Args:
        material (Material): Material object
        atom_index (int): Index of the central atom
        coordinate (List[float]): Coordinate of the central atom in crystal coordinates
        layer_thickness (float): Thickness of the layer in angstroms

    Returns:
        List[int]: List of indices of atoms within the specified layer
    """
    if coordinate is not None:
        return get_atom_indices_within_layer_by_atom_coordinate(material, coordinate, layer_thickness)
    if atom_index is not None:
        return get_atom_indices_within_layer_by_atom_index(material, atom_index, layer_thickness)


def get_atom_indices_within_radius_pbc(
    material: Material, atom_index: Optional[int] = 0, coordinate: Optional[List[float]] = None, radius: float = 1
):
    """
    Select all atoms within a specified radius of a central atom considering periodic boundary conditions.

    Args:
        material (Material): Material object
        atom_index (int): Index of the central atom
        coordinate (List[float]): Coordinate of the central atom in crystal coordinates
        radius (float): Radius of the sphere in angstroms

    Returns:
        List[int]: List of indices of atoms within the specified
    """

    if coordinate is not None:
        atom_index = get_closest_site_id_from_coordinate(material, coordinate)

    structure = to_pymatgen(material)
    immutable_structure = PymatgenIStructure.from_sites(structure.sites)

    central_atom = immutable_structure[atom_index]
    sites_within_radius = structure.get_sites_in_sphere(central_atom.coords, radius)

    selected_indices = [site.index for site in sites_within_radius]
    return selected_indices


def get_atom_indices_with_condition_on_coordinates(
    material: Material,
    condition: Callable[[List[float]], bool],
    use_cartesian_coordinates: bool = False,
) -> List[int]:
    """
    Select atoms whose coordinates satisfy the given condition.

    Args:
        material (Material): Material object
        condition (Callable[List[float], bool]): Function that checks if coordinates satisfy the condition.
        use_cartesian_coordinates (bool): Whether to use Cartesian coordinates for the condition evaluation.

    Returns:
        List[int]: List of indices of atoms whose coordinates satisfy the condition.
    """
    new_material = material.clone()
    new_basis = new_material.basis
    if use_cartesian_coordinates:
        new_basis.to_cartesian()
    else:
        new_basis.to_crystal()
    new_material.basis = new_basis
    coordinates = new_material.basis.coordinates.to_array_of_values_with_ids()

    selected_indices = []
    for coord in coordinates:
        if condition(coord.value):
            selected_indices.append(coord.id)

    return selected_indices


def get_nearest_neighbors_atom_indices(
    material: Material,
    coordinate: Optional[List[float]] = None,
    cutoff: float = 15.0,
) -> Optional[List[int]]:
    """
    Returns the indices of direct neighboring atoms to a specified position in the material using Voronoi tessellation.

    Args:
        material (Material): The material object to find neighbors in.
        coordinate (List[float]): The position to find neighbors for.
        cutoff (float): The cutoff radius for identifying neighbors.

    Returns:
        List[int]: A list of indices of neighboring atoms, or an empty list if no neighbors are found.
    """
    if coordinate is None:
        coordinate = [0, 0, 0]
    structure = to_pymatgen(material)
    voronoi_nn = PymatgenVoronoiNN(
        tol=0.1,
        cutoff=cutoff,
        allow_pathological=False,
        weight="solid_angle",
        extra_nn_info=True,
        compute_adj_neighbors=False,
    )
    coordinates = material.basis.coordinates
    coordinates.filter_by_values([coordinate])
    site_index = coordinates.ids[0]
    remove_added = False
    if site_index is None:
        structure.append("X", coordinate, validate_proximity=False)
        site_index = len(structure.sites) - 1
        remove_added = True
    neighbors = voronoi_nn.get_nn_info(structure, site_index)
    neighboring_atoms_pymatgen_ids = [n["site_index"] for n in neighbors]

    if remove_added:
        structure.remove_sites([-1])

    all_coordinates = material.basis.coordinates
    all_coordinates.filter_by_indices(neighboring_atoms_pymatgen_ids)
    return all_coordinates.ids


def get_atomic_coordinates_extremum(
    material: Material,
    extremum: Literal["max", "min"] = "max",
    axis: Literal["x", "y", "z"] = "z",
    use_cartesian_coordinates: bool = False,
) -> float:
    """
    Return minimum or maximum of coordinates along the specified axis.

    Args:
        material (Material): Material object.
        extremum (str): "min" or "max".
        axis (str):  "x", "y", or "z".
        use_cartesian_coordinates (bool): Whether to use Cartesian coordinates.
    Returns:
        float: Minimum or maximum of coordinates along the specified axis.
    """
    new_material = material.clone()
    new_basis = new_material.basis
    if use_cartesian_coordinates:
        new_basis.to_cartesian()
    else:
        new_basis.to_crystal()
    new_material.basis = new_basis
    coordinates = new_material.basis.coordinates.to_array_of_values_with_ids()
    values = [coord.value[{"x": 0, "y": 1, "z": 2}[axis]] for coord in coordinates]
    return getattr(np, extremum)(values)


def get_undercoordinated_atom_indices(
    material: Material, indices_to_check: List[int]
) -> Tuple[List[int], List[List[int]]]:
    neighbors_indices_array = []
    neighbors_numbers = []
    undercoordinated_atom_indices: List[int] = []
    set_of_neighbors_numbers = set()
    for idx in indices_to_check:
        coordinate = material.basis.coordinates.values[idx]
        neighbors_indices: List[int] = []
        try:
            neighbors_indices = get_nearest_neighbors_atom_indices(material, coordinate, cutoff=5)
            neighbors_indices_array.append(neighbors_indices)
        except Exception as e:
            print(f"Error: {e}")
            neighbors_indices_array.append([])
            continue
        neighbors_numbers.append(len(neighbors_indices))
        set_of_neighbors_numbers.add(len(neighbors_indices))
    threshold = max(set_of_neighbors_numbers)
    for idx, number in enumerate(neighbors_numbers):
        if number < threshold:
            undercoordinated_atom_indices.append(idx)
    return undercoordinated_atom_indices, neighbors_indices_array
