import numpy as np
from ase import Atoms

from .convert import decorator_convert_material_args_kwargs_to_atoms


@decorator_convert_material_args_kwargs_to_atoms
def get_average_interlayer_distance(
    interface_atoms: Atoms, tag_substrate: str, tag_film: str, threshold: float = 0.5
) -> float:
    """
    Calculate the average distance between the top layer of substrate atoms and the bottom layer of film atoms.

    Args:
        interface_atoms (ase.Atoms): The ASE Atoms object containing both sets of atoms.
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
def get_surface_area(atoms: Atoms):
    """
    Calculate the area of the surface perpendicular to the z-axis of the atoms structure.

    Args:
        atoms (ase.Atoms): The Atoms object to calculate the surface area of.

    Returns:
        float: The surface area of the atoms.
    """
    matrix = atoms.cell
    cross_product = np.cross(matrix[0], matrix[1])
    return np.linalg.norm(cross_product)
