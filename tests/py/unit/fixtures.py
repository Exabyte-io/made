from typing import Tuple

from ase.build import bulk
from mat3ra.made.material import Material
from mat3ra.made.tools.build.interface import interface_patch_with_mean_abs_strain
from mat3ra.made.tools.convert import from_ase
from pymatgen.analysis.elasticity.strain import Strain
from pymatgen.core.interface import Interface

from .utils import atoms_to_interface_structure

# ASE Atoms fixtures
substrate = bulk("Si", cubic=True)
film = bulk("Cu", cubic=True)
INTERFACE_ATOMS = substrate + film
INTERFACE_ATOMS.set_tags([1] * len(substrate) + [2] * len(film))

# Material fixtures
SUBSTRATE_MATERIAL = Material(from_ase(substrate))
LAYER_MATERIAL = Material(from_ase(film))

# Pymatgen Interface fixtures
INTERFACE_TERMINATION: Tuple = ("Si_termination", "Cu_termination")

interface_structure = atoms_to_interface_structure(INTERFACE_ATOMS)
dict = interface_structure.as_dict()

INTERFACE_STRUCTURE = Interface.from_dict(dict)
# Create properties that are assigned during interface creation in ZSL algorithm as dict and json
INTERFACE_PROPERTIES_MOCK = {
    "film_sl_vectors": [[2.0, -4.0, 0.0], [-2.0, -4.0, 0.0]],
    "substrate_sl_vectors": [[-3.5, 3.5, 0.0], [-3.5, 0.0, 3.5]],
    "film_vectors": [[1.234, -2.345, 0.0], [1.234, 2.345, 0.0]],
    "substrate_vectors": [[-3.5, 3.5, 0.0], [-3.5, 0.0, 3.5]],
    "film_transformation": [[2.0, 0.0], [0.0, 2.0]],
    "substrate_transformation": [[1.0, 0.0], [0.0, 1.0]],
    "strain": Strain(
        [
            [4.74636368e-03, -1.07111363e-16, -2.13208193e-17],
            [-1.07111363e-16, 4.74636368e-03, 5.36172802e-18],
            [-2.13208193e-17, 5.36172802e-18, -5.55111512e-17],
        ]
    ),
    "von_mises_strain": 0.1,
    "termination": INTERFACE_TERMINATION,
    "film_thickness": 1,
    "substrate_thickness": 3,
}
INTERFACE_PROPERTIES_JSON = {
    "film_sl_vectors": [[2.0, -4.0, 0.0], [-2.0, -4.0, 0.0]],
    "substrate_sl_vectors": [[-3.5, 3.5, 0.0], [-3.5, 0.0, 3.5]],
    "film_vectors": [[1.234, -2.345, 0.0], [1.234, 2.345, 0.0]],
    "substrate_vectors": [[-3.5, 3.5, 0.0], [-3.5, 0.0, 3.5]],
    "film_transformation": [[2.0, 0.0], [0.0, 2.0]],
    "substrate_transformation": [[1.0, 0.0], [0.0, 1.0]],
    "strain": [
        [0.0047463636820092026, -1.0711136330719372e-16, -2.132081932243309e-17],
        [-1.0711136330719372e-16, 0.0047463636820092026, 5.3617280208039476e-18],
        [-2.132081932243309e-17, 5.3617280208039476e-18, -5.551115123125783e-17],
    ],
    "von_mises_strain": 0.1,
    "termination": INTERFACE_TERMINATION,
    "film_thickness": 1,
    "substrate_thickness": 3,
    "mean_abs_strain": 0.0010500000000000002,
}

# Add properties to interface structure
INTERFACE_STRUCTURE.interface_properties = INTERFACE_PROPERTIES_MOCK
INTERFACE_STRUCTURE = interface_patch_with_mean_abs_strain(INTERFACE_STRUCTURE)
