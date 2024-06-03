import inspect
import json
from functools import wraps
from typing import Any, Callable, Dict, Union

from ase import Atoms
from mat3ra.utils.mixins import RoundNumericValuesMixin
from mat3ra.utils.object import NumpyNDArrayRoundEncoder
from pymatgen.core.interface import Interface, label_termination
from pymatgen.core.structure import Lattice, Structure
from pymatgen.core.surface import Slab
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.vasp.inputs import Poscar

from ..material import Material

PymatgenStructure = Structure
PymatgenSlab = Slab
PymatgenInterface = Interface
ASEAtoms = Atoms
label_pymatgen_slab_termination = label_termination


def to_pymatgen(material_or_material_data: Union[Material, Dict[str, Any]]) -> Structure:
    """
    Converts material object in ESSE format to a pymatgen Structure object.

    Args:
        material_or_material_data (dict|class): A dictionary containing the material information in ESSE format.

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
    lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)

    basis = material_data["basis"]
    elements = [element["value"] for element in basis["elements"]]
    coordinates = [coord["value"] for coord in basis["coordinates"]]
    labels = [label["value"] for label in basis.get("labels", [])]

    # Assuming that the basis units are fractional since it's a crystal basis
    coords_are_cartesian = "units" in basis and basis["units"].lower() == "angstrom"

    if "labels" in inspect.signature(Structure.__init__).parameters:
        structure = Structure(lattice, elements, coordinates, coords_are_cartesian=coords_are_cartesian, labels=labels)
    else:
        # Passing labels does not work for pymatgen `2023.6.23` supporting py3.8
        print(f"labels: {labels}. Not passing labels to pymatgen.")
        structure = Structure(lattice, elements, coordinates, coords_are_cartesian=coords_are_cartesian)

    return structure


def from_pymatgen(structure: Union[Structure, Interface]) -> Dict[str, Any]:
    """
    Converts a pymatgen Structure object to a material object in ESSE format.

    Args:
        structure (Structure): A pymatgen Structure object.

    Returns:
        dict: A dictionary containing the material information in ESSE format.
    """
    __round__ = RoundNumericValuesMixin.round_array_or_number
    basis = {
        "elements": [{"id": i, "value": str(site.specie)} for i, site in enumerate(structure.sites)],
        "coordinates": [
            {"id": i, "value": __round__(list(site.frac_coords))} for i, site in enumerate(structure.sites)
        ],
        "units": "crystal",
        "cell": __round__(structure.lattice.matrix.tolist()),
        "constraints": [],
    }

    lattice = {
        "a": __round__(structure.lattice.a),
        "b": __round__(structure.lattice.b),
        "c": __round__(structure.lattice.c),
        "alpha": __round__(structure.lattice.alpha),
        "beta": __round__(structure.lattice.beta),
        "gamma": __round__(structure.lattice.gamma),
        "units": {"length": "angstrom", "angle": "degree"},
        "type": "TRI",
        "vectors": {
            "a": __round__(structure.lattice.matrix[0].tolist()),
            "b": __round__(structure.lattice.matrix[1].tolist()),
            "c": __round__(structure.lattice.matrix[2].tolist()),
            "alat": 1,
            "units": "angstrom",
        },
    }

    metadata = {"boundaryConditions": {"type": "pbc", "offset": 0}}

    # TODO: consider using Interface JSONSchema from ESSE when such created and adapt interface_properties accordingly.
    # Add interface properties to metadata according to pymatgen Interface as a JSON object
    if hasattr(structure, "interface_properties"):
        interface_props = structure.interface_properties
        # TODO: figure out how to round the values and stringify terminations tuple
        #  in the interface properties with Encoder
        for key, value in interface_props.items():
            if isinstance(value, tuple):
                interface_props[key] = str(value)
        metadata["interface_properties"] = json.loads(json.dumps(interface_props, cls=NumpyNDArrayRoundEncoder))

    material_data = {
        "name": structure.formula,
        "basis": basis,
        "lattice": lattice,
        "isNonPeriodic": not structure.is_ordered,
        "_id": "",
        "metadata": metadata,
        "isUpdated": True,
    }

    return material_data


def to_poscar(material_or_material_data: Union[Material, Dict[str, Any]]) -> str:
    """
    Converts material object in ESSE format to a POSCAR string.

    Args:
        material_or_material_data (dict|class): A dictionary containing the material information in ESSE format.

    Returns:
        str: A POSCAR string.
    """
    structure = to_pymatgen(material_or_material_data)
    poscar = Poscar(structure)
    # For pymatgen `2023.6.23` supporting py3.8 the method name is "get_string"
    # TODO: cleanup the if statement when dropping support for py3.8
    if hasattr(poscar, "get_string"):
        return poscar.get_string()
    return poscar.get_str()


def from_poscar(poscar: str) -> Dict[str, Any]:
    """
    Converts a POSCAR string to a material object in ESSE format.

    Args:
        poscar (str): A POSCAR string.

    Returns:
        dict: A dictionary containing the material information in ESSE format.
    """
    structure = Structure.from_str(poscar, "poscar")
    return from_pymatgen(structure)


def to_ase(material_or_material_data: Union[Material, Dict[str, Any]]) -> Atoms:
    """
    Converts material object in ESSE format to an ASE Atoms object.

    Args:
        material_or_material_data (dict|class): A dictionary containing the material information in ESSE format.

    Returns:
        Any: An ASE Atoms object.
    """
    # TODO: check that atomic labels are properly handled
    structure = to_pymatgen(material_or_material_data)
    return AseAtomsAdaptor.get_atoms(structure)


def from_ase(ase_atoms: Atoms) -> Dict[str, Any]:
    """
    Converts an ASE Atoms object to a material object in ESSE format.

    Args:
        ase_atoms (Atoms): An ASE Atoms object.

    Returns:
        dict: A dictionary containing the material information in ESSE format.
    """
    # TODO: check that atomic labels/tags are properly handled
    structure = AseAtomsAdaptor.get_structure(ase_atoms)
    return from_pymatgen(structure)


def decorator_convert_material_args_kwargs_to_atoms(func: Callable) -> Callable:
    """
    Decorator that converts ESSE Material objects to ASE Atoms objects.
    """

    @wraps(func)
    def wrapper(*args, **kwargs):
        # Convert args if they are of type ESSE Material
        new_args = [to_ase(arg) if isinstance(arg, Material) else arg for arg in args]

        # Convert kwargs if they are of type ESSE Material
        new_kwargs = {k: to_ase(v) if isinstance(v, Material) else v for k, v in kwargs.items()}

        # Call the original function with the converted arguments
        return func(*new_args, **new_kwargs)

    return wrapper


def decorator_convert_material_args_kwargs_to_structure(func: Callable) -> Callable:
    """
    Decorator that converts ESSE Material objects to pymatgen Structure objects.
    """

    @wraps(func)
    def wrapper(*args, **kwargs):
        # Convert args if they are of type ESSE Material
        new_args = [to_pymatgen(arg) if isinstance(arg, Material) else arg for arg in args]

        # Convert kwargs if they are of type ESSE Material
        new_kwargs = {k: to_pymatgen(v) if isinstance(v, Material) else v for k, v in kwargs.items()}

        # Call the original function with the converted arguments
        return func(*new_args, **new_kwargs)

    return wrapper


def convert_atoms_or_structure_to_material(item):
    if isinstance(item, Structure):
        return from_pymatgen(item)
    elif isinstance(item, Atoms):
        return from_ase(item)
    return item
