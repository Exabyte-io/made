from typing import Callable, List, Literal, Optional

import numpy as np
from mat3ra.made.material import Material
from scipy.spatial import cKDTree

from .convert import decorator_convert_material_args_kwargs_to_atoms, to_pymatgen
from .enums import SurfaceTypes
from .third_party import ASEAtoms, PymatgenIStructure, PymatgenVoronoiNN
from .utils import decorator_handle_periodic_boundary_conditions


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


def is_height_within_limits(z: float, z_extremum: float, depth: float, surface: SurfaceTypes) -> bool:
    """
    Check if the height of an atom is within the specified limits.

    Args:
        z (float): The z-coordinate of the atom.
        z_extremum (float): The extremum z-coordinate of the surface.
        depth (float): The depth from the surface to look for exposed atoms.
        surface (SurfaceTypes): The surface type (top or bottom).

    Returns:
        bool: True if the height is within the limits, False otherwise.
    """
    return (z >= z_extremum - depth) if surface == SurfaceTypes.TOP else (z <= z_extremum + depth)


def is_shadowed_by_neighbors_from_surface(
    z: float, neighbors_indices: List[int], surface: SurfaceTypes, coordinates: np.ndarray
) -> bool:
    """
    Check if any one of the neighboring atoms shadow the atom from the surface by being closer to the specified surface.

    Args:
        z (float): The z-coordinate of the atom.
        neighbors_indices (List[int]): List of indices of neighboring atoms.
        surface (SurfaceTypes): The surface type (top or bottom).
        coordinates (np.ndarray): The coordinates of the atoms.

    Returns:
        bool: True if the atom is not shadowed, False otherwise.
    """
    return not any(
        (coordinates[n][2] > z if surface == SurfaceTypes.TOP else coordinates[n][2] < z) for n in neighbors_indices
    )


@decorator_handle_periodic_boundary_conditions(cutoff=0.1)
def get_surface_atom_indices(
    material: Material, surface: SurfaceTypes = SurfaceTypes.TOP, shadowing_radius: float = 2.5, depth: float = 5
) -> List[int]:
    """
    Identify exposed atoms on the top or bottom surface of the material.

    Args:
        material (Material): Material object to get surface atoms from.
        surface (SurfaceTypes): Specify "top" or "bottom" to detect the respective surface atoms.
        shadowing_radius (float): Radius for atoms shadowing underlying from detecting as exposed.
        depth (float): Depth from the surface to look for exposed atoms.

    Returns:
        List[int]: List of indices of exposed surface atoms.
    """
    new_material = material.clone()
    new_material.to_cartesian()
    coordinates = np.array(new_material.basis.coordinates.values)
    ids = new_material.basis.coordinates.ids
    kd_tree = cKDTree(coordinates)

    z_extremum = np.max(coordinates[:, 2]) if surface == SurfaceTypes.TOP else np.min(coordinates[:, 2])

    exposed_atoms_indices = []
    for idx, (x, y, z) in enumerate(coordinates):
        if is_height_within_limits(z, z_extremum, depth, surface):
            neighbors_indices = kd_tree.query_ball_point([x, y, z], r=shadowing_radius)
            if is_shadowed_by_neighbors_from_surface(z, neighbors_indices, surface, coordinates):
                exposed_atoms_indices.append(ids[idx])

    return exposed_atoms_indices


def get_coordination_numbers(
    material: Material,
    indices: Optional[List[int]] = None,
    cutoff: float = 3.0,
) -> List[int]:
    """
    Calculate the coordination numbers of atoms in the material.

    Args:
        material (Material): Material object to calculate coordination numbers for.
        indices (List[int]): List of atom indices to calculate coordination numbers for.
        cutoff (float): The cutoff radius for identifying neighbors.

    Returns:
        List[int]: List of coordination numbers for each atom in the material.
    """
    new_material = material.clone()
    new_material.to_cartesian()
    if indices is not None:
        new_material.basis.coordinates.filter_by_indices(indices)
    coordinates = np.array(new_material.basis.coordinates.values)
    kd_tree = cKDTree(coordinates)

    coordination_numbers = []
    for idx, (x, y, z) in enumerate(coordinates):
        neighbors = kd_tree.query_ball_point([x, y, z], r=cutoff)
        # Explicitly remove the atom itself from the list of neighbors
        neighbors = [n for n in neighbors if n != idx]
        coordination_numbers.append(len(neighbors))

    return coordination_numbers


@decorator_handle_periodic_boundary_conditions(cutoff=0.1)
def get_undercoordinated_atom_indices(
    material: Material,
    indices: List[int],
    cutoff: float = 3.0,
    coordination_threshold: int = 3,
) -> List[int]:
    """
    Identify undercoordinated atoms among the specified indices in the material.

    Args:
        material (Material): Material object to identify undercoordinated atoms in.
        indices (List[int]): List of atom indices to check for undercoordination.
        cutoff (float): The cutoff radius for identifying neighbors.
        coordination_threshold (int): The coordination number threshold for undercoordination.

    Returns:
        List[int]: List of indices of undercoordinated atoms.
    """
    coordination_numbers = get_coordination_numbers(material, indices, cutoff)
    undercoordinated_atoms_indices = [i for i, cn in enumerate(coordination_numbers) if cn <= coordination_threshold]
    return undercoordinated_atoms_indices


def get_local_extremum_atom_index(
    material: Material,
    coordinate: List[float],
    extremum: Literal["max", "min"] = "max",
    vicinity: float = 1.0,
    use_cartesian_coordinates: bool = False,
) -> int:
    """
    Return the id of the atom with the minimum or maximum z-coordinate
    within a certain vicinity of a given (x, y) coordinate.

    Args:
        material (Material): Material object.
        coordinate (List[float]): (x, y, z) coordinate to find the local extremum atom index for.
        extremum (str): "min" or "max".
        vicinity (float): Radius of the vicinity, in Angstroms.
        use_cartesian_coordinates (bool): Whether to use Cartesian coordinates.

    Returns:
        int: id of the atom with the minimum or maximum z-coordinate.
    """
    new_material = material.clone()
    new_material.to_cartesian()
    if not use_cartesian_coordinates:
        coordinate = new_material.basis.cell.convert_point_to_cartesian(coordinate)

    coordinates = np.array(new_material.basis.coordinates.values)
    ids = np.array(new_material.basis.coordinates.ids)
    tree = cKDTree(coordinates[:, :2])
    indices = tree.query_ball_point(coordinate[:2], vicinity)
    z_values = [(id, coord[2]) for id, coord in zip(ids[indices], coordinates[indices])]

    if extremum == "max":
        extremum_z_atom = max(z_values, key=lambda item: item[1])
    else:
        extremum_z_atom = min(z_values, key=lambda item: item[1])

    return extremum_z_atom[0]
