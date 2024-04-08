from typing import Any, Dict, Union

from mat3ra.made.material import Material
from pymatgen.core.structure import Lattice, Structure


def to_pymatgen(material_or_material_data: Union[Material, Dict[str, Any]]) -> Structure:
    """
    Converts material object in ESSE format to a pymatgen Structure object.

    Args:
        material_data (dict): A dictionary containing the material information in ESSE format.

    Returns:
        Structure: A pymatgen Structure object.
    """
    material_data = material_or_material_data

    if isinstance(material_or_material_data, Material):
        material_data = material_or_material_data.to_json()

    lattice_params = material_data["lattice"]
    a = lattice_params["a"]
    b = lattice_params["b"]
    c = lattice_params["c"]
    alpha = lattice_params["alpha"]
    beta = lattice_params["beta"]
    gamma = lattice_params["gamma"]

    # Create a Lattice from parameters
    lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)

    # Extract the basis information
    basis = material_data["basis"]
    elements = [element["value"] for element in basis["elements"]]
    coordinates = [coord["value"] for coord in basis["coordinates"]]

    # Assuming that the basis units are fractional since it's a crystal basis
    coords_are_cartesian = "units" in basis and basis["units"].lower() == "angstrom"

    # Create the Structure
    structure = Structure(lattice, elements, coordinates, coords_are_cartesian=coords_are_cartesian)

    return structure
