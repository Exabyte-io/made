from typing import List, Optional

import numpy as np

from ..material import Material
from .convert import decorator_convert_material_args_kwargs_to_atoms, to_pymatgen
from .third_party import ASEAtoms, PymatgenIStructure


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


def get_closest_site_id_from_position(material: Material, position: List[float]) -> int:
    """
    Get the site ID of the closest site to a given position in the crystal.

    Args:
        material (Material): The material object to find the closest site in.
        position (List[float]): The position to find the closest site to.

    Returns:
        int: The site ID of the closest site.
    """
    coordinates = np.array(material.coordinates_array)
    position = np.array(position)  # type: ignore
    distances = np.linalg.norm(coordinates - position, axis=1)
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
    central_atom_position = coordinates[atom_index]
    central_atom_projection = np.dot(central_atom_position.value, direction_norm)

    layer_thickness_frac = layer_thickness / direction_length

    lower_bound = central_atom_projection - layer_thickness_frac / 2
    upper_bound = central_atom_projection + layer_thickness_frac / 2

    selected_indices = []
    for coord in coordinates:
        # Project each position onto the direction vector
        projection = np.dot(coord.value, direction_norm)
        if lower_bound <= projection <= upper_bound:
            selected_indices.append(coord.id)
    return selected_indices


def get_atom_indices_within_layer_by_atom_position(material: Material, position: List[float], layer_thickness: float):
    """
    Select all atoms within a specified layer thickness of a central atom along a direction.
    This direction will be orthogonal to the AB plane.
    Layer thickness is converted from angstroms to fractional units based on the lattice vector length.

    Args:
        material (Material): Material object
        position (List[float]): Position of the central atom in crystal coordinates
        layer_thickness (float): Thickness of the layer in angstroms

    Returns:
        List[int]: List of indices of atoms within the specified layer
    """
    site_id = get_closest_site_id_from_position(material, position)
    return get_atom_indices_within_layer_by_atom_index(material, site_id, layer_thickness)


def get_atom_indices_within_layer(
    material: Material,
    atom_index: Optional[int] = 0,
    position: Optional[List[float]] = None,
    layer_thickness: float = 1,
):
    """
    Select all atoms within a specified layer thickness of the central atom along the c-vector direction.

    Args:
        material (Material): Material object
        atom_index (int): Index of the central atom
        position (List[float]): Position of the central atom in crystal coordinates
        layer_thickness (float): Thickness of the layer in angstroms

    Returns:
        List[int]: List of indices of atoms within the specified layer
    """
    if position is not None:
        return get_atom_indices_within_layer_by_atom_position(material, position, layer_thickness)
    if atom_index is not None:
        return get_atom_indices_within_layer_by_atom_index(material, atom_index, layer_thickness)


def get_atom_indices_within_radius_pbc(
    material: Material, atom_index: Optional[int] = 0, position: Optional[List[float]] = None, radius: float = 1
):
    """
    Select all atoms within a specified radius of a central atom considering periodic boundary conditions.

    Args:
        material (Material): Material object
        atom_index (int): Index of the central atom
        position (List[float]): Position of the central atom in crystal coordinates
        radius (float): Radius of the sphere in angstroms

    Returns:
        List[int]: List of indices of atoms within the specified
    """

    if position is not None:
        atom_index = get_closest_site_id_from_position(material, position)

    structure = to_pymatgen(material)
    immutable_structure = PymatgenIStructure.from_sites(structure.sites)

    central_atom = immutable_structure[atom_index]
    sites_within_radius = structure.get_sites_in_sphere(central_atom.coords, radius)

    selected_indices = [site.index for site in sites_within_radius]
    return selected_indices
