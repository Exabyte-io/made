from ase.build import bulk
from mat3ra.made.material import Material
from mat3ra.made.tools.build.interface.termination_pair import TerminationPair
from mat3ra.made.tools.build.slab import SlabConfiguration, get_terminations
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

SUBSTRATE_CONFIGURATION = SlabConfiguration(bulk=SUBSTRATE_MATERIAL, thickness=3)
LAYER_CONFIGURATION = SlabConfiguration(bulk=LAYER_MATERIAL)

layer_terminations = get_terminations(LAYER_CONFIGURATION)
substrate_terminations = get_terminations(SUBSTRATE_CONFIGURATION)

# Pymatgen Interface fixtures
INTERFACE_TERMINATION_PAIR: TerminationPair = TerminationPair(
    (
        layer_terminations[0],
        substrate_terminations[0],
    )
)
INTERFACE_TERMINATION_AS_STR = str(INTERFACE_TERMINATION_PAIR.self)

interface_structure = atoms_to_interface_structure(INTERFACE_ATOMS)
dict = interface_structure.as_dict()

INTERFACE_STRUCTURE = Interface.from_dict(dict)
# Create properties that are assigned during interface creation in ZSL algorithm as dict and json
INTERFACE_PROPERTIES_MOCK = {
    "film_transformation": [[2.0, 0.0], [0.0, 2.0]],
    "substrate_transformation": [[1.0, 0.0], [0.0, 1.0]],
    "strain": Strain([[0.004746364, -0.0, -0.0], [-0.0, 0.004746364, 0.0], [-0.0, 0.0, -0.0]]),
    "von_mises_strain": 0.001,
    "termination": (INTERFACE_TERMINATION_PAIR.film_termination, INTERFACE_TERMINATION_PAIR.substrate_termination),
}
INTERFACE_PROPERTIES_JSON = {
    "film_transformation": [[2.0, 0.0], [0.0, 2.0]],
    "substrate_transformation": [[1.0, 0.0], [0.0, 1.0]],
    "strain": [[0.004746364, -0.0, -0.0], [-0.0, 0.004746364, 0.0], [-0.0, 0.0, -0.0]],
    "von_mises_strain": 0.001,
    "termination": INTERFACE_TERMINATION_AS_STR,
    # "mean_abs_strain": 0.00105,
}


# Add properties to interface structure
INTERFACE_STRUCTURE.interface_properties = INTERFACE_PROPERTIES_MOCK

# TODO: Use fixtures package when available
SI_SUPERCELL_2X2X1 = {
    "name": "Si8",
    "basis": {
        "elements": [
            {"id": 0, "value": "Si"},
            {"id": 1, "value": "Si"},
            {"id": 2, "value": "Si"},
            {"id": 3, "value": "Si"},
            {"id": 4, "value": "Si"},
            {"id": 5, "value": "Si"},
            {"id": 6, "value": "Si"},
            {"id": 7, "value": "Si"},
        ],
        "coordinates": [
            {"id": 0, "value": [0.0, 0.0, 0.0]},
            {"id": 1, "value": [0.125, 0.125, 0.25]},
            {"id": 2, "value": [0.0, 0.5, 0.0]},
            {"id": 3, "value": [0.125, 0.625, 0.25]},
            {"id": 4, "value": [0.5, 0.0, 0.0]},
            {"id": 5, "value": [0.625, 0.125, 0.25]},
            {"id": 6, "value": [0.5, 0.5, 0.0]},
            {"id": 7, "value": [0.625, 0.625, 0.25]},
        ],
        "units": "crystal",
        "cell": [[6.697840473, 0.0, 3.867], [2.232613491, 6.314784557, 3.867], [0.0, 0.0, 3.867]],
        "constraints": [],
    },
    "lattice": {
        "a": 7.734,
        "b": 7.734,
        "c": 3.867,
        "alpha": 60.0,
        "beta": 60.0,
        "gamma": 60.0,
        "units": {"length": "angstrom", "angle": "degree"},
        "type": "TRI",
        "vectors": {
            "a": [6.697840473, 0.0, 3.867],
            "b": [2.232613491, 6.314784557, 3.867],
            "c": [0.0, 0.0, 3.867],
            "alat": 1,
            "units": "angstrom",
        },
    },
    "isNonPeriodic": False,
    "_id": "",
    "metadata": {"boundaryConditions": {"type": "pbc", "offset": 0}},
    "isUpdated": True,
}

SI_SLAB = {
    "name": "Si8(001), termination Si_P4/mmm_1, Slab",
    "basis": {
        "elements": [
            {"id": 0, "value": "Si"},
            {"id": 1, "value": "Si"},
            {"id": 2, "value": "Si"},
            {"id": 3, "value": "Si"},
        ],
        "coordinates": [
            {"id": 0, "value": [0.5, 0.5, 0.5625]},
            {"id": 1, "value": [0.5, 0.0, 0.6875]},
            {"id": 2, "value": [0.0, 0.0, 0.8125]},
            {"id": 3, "value": [-0.0, 0.5, 0.9375]},
        ],
        "units": "crystal",
        "cell": [[3.867, 0.0, 0.0], [0.0, 3.867, 0.0], [0.0, 0.0, 10.937527692]],
        "constraints": [],
    },
    "lattice": {
        "a": 3.867,
        "b": 3.867,
        "c": 10.937527692,
        "alpha": 90.0,
        "beta": 90.0,
        "gamma": 90.0,
        "units": {"length": "angstrom", "angle": "degree"},
        "type": "TRI",
        "vectors": {
            "a": [3.867, 0.0, 0.0],
            "b": [0.0, 3.867, 0.0],
            "c": [0.0, 0.0, 10.937527692],
            "alat": 1,
            "units": "angstrom",
        },
    },
    "isNonPeriodic": False,
    "_id": "",
    "metadata": {
        "boundaryConditions": {"type": "pbc", "offset": 0},
        "termination": "Si_P4/mmm_1",
    },
    "isUpdated": True,
}
