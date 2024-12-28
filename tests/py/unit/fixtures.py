import copy
from functools import reduce
from typing import Any, Dict

from ase.build import bulk
from mat3ra.made.material import Material
from mat3ra.made.tools.build.interface.termination_pair import TerminationPair
from mat3ra.made.tools.build.slab import SlabConfiguration, create_slab, get_terminations
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
FILM_MATERIAL = Material(from_ase(film))

SUBSTRATE_CONFIGURATION = SlabConfiguration(bulk=SUBSTRATE_MATERIAL, thickness=3)
FILM_CONFIGURATION = SlabConfiguration(bulk=FILM_MATERIAL)

substrate_terminations = get_terminations(SUBSTRATE_CONFIGURATION)
film_terminations = get_terminations(FILM_CONFIGURATION)

# Pymatgen Interface fixtures
INTERFACE_TERMINATION_PAIR: TerminationPair = TerminationPair(
    film_terminations[0],
    substrate_terminations[0],
)
INTERFACE_TERMINATION_AS_STR = str(INTERFACE_TERMINATION_PAIR)

interface_structure = atoms_to_interface_structure(INTERFACE_ATOMS)
interface_dict = interface_structure.as_dict()

INTERFACE_STRUCTURE = Interface.from_dict(interface_dict)
# Create properties that are assigned during interface creation in ZSL algorithm as dict and json
INTERFACE_PROPERTIES_MOCK = {
    "film_transformation": [[2.0, 0.0], [0.0, 2.0]],
    "substrate_transformation": [[1.0, 0.0], [0.0, 1.0]],
    "strain": Strain([[0.004746364, -0.0, -0.0], [-0.0, 0.004746364, 0.0], [-0.0, 0.0, -0.0]]),
    "von_mises_strain": 0.001,
    "mean_abs_strain": 0.00105,
    "termination": (INTERFACE_TERMINATION_PAIR.film_termination, INTERFACE_TERMINATION_PAIR.substrate_termination),
}
INTERFACE_PROPERTIES_JSON = {
    "film_transformation": [[2.0, 0.0], [0.0, 2.0]],
    "substrate_transformation": [[1.0, 0.0], [0.0, 1.0]],
    "strain": [[0.004746364, -0.0, -0.0], [-0.0, 0.004746364, 0.0], [-0.0, 0.0, -0.0]],
    "von_mises_strain": 0.001,
    "termination": INTERFACE_TERMINATION_AS_STR,
    "mean_abs_strain": 0.00105,
}

# Add properties to interface structure
INTERFACE_STRUCTURE.interface_properties = INTERFACE_PROPERTIES_MOCK
INTERFACE_NAME = "Cu4(001)-Si8(001), Interface, Strain 0.062pct"

# TODO: Use fixtures package when available
SI_CONVENTIONAL_CELL: Dict[str, Any] = {
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
            {"id": 0, "value": [0.5, 0.0, 0.0]},
            {"id": 1, "value": [0.25, 0.25, 0.75]},
            {"id": 2, "value": [0.5, 0.5, 0.5]},
            {"id": 3, "value": [0.25, 0.75, 0.25]},
            {"id": 4, "value": [0.0, 0.0, 0.5]},
            {"id": 5, "value": [0.75, 0.25, 0.25]},
            {"id": 6, "value": [0.0, 0.5, 0.0]},
            {"id": 7, "value": [0.75, 0.75, 0.75]},
        ],
        "units": "crystal",
        "constraints": [],
        "labels": [],
    },
    "lattice": {
        "a": 5.468763846,
        "b": 5.468763846,
        "c": 5.468763846,
        "alpha": 90.0,
        "beta": 90.0,
        "gamma": 90.0,
        "units": {"length": "angstrom", "angle": "degree"},
        "type": "TRI",
        "vectors": {
            "a": [5.468763846, 0.0, 0.0],
            "b": [-0.0, 5.468763846, 0.0],
            "c": [0.0, 0.0, 5.468763846],
            "alat": 1,
            "units": "angstrom",
        },
    },
    "isNonPeriodic": False,
    "_id": "",
    "metadata": {"boundaryConditions": {"type": "pbc", "offset": 0}},
    "isUpdated": True,
}

SI_SUPERCELL_2X2X1: Dict[str, Any] = {
    "name": "Silicon FCC",
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
        "constraints": [],
        "labels": [],
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

SI_SLAB_CONFIGURATION: Dict[str, Any] = {
    "type": "SlabConfiguration",
    "bulk": SI_CONVENTIONAL_CELL,
    "miller_indices": (0, 0, 1),
    "thickness": 1,
    "vacuum": 5.0,
    "xy_supercell_matrix": [[1, 0], [0, 1]],
    "use_conventional_cell": True,
    "use_orthogonal_z": True,
    "make_primitive": True,
}

SI_SLAB_100: Dict[str, Any] = {
    "name": "Si8(001), termination Si_P4/mmm_1, Slab",
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
            {"id": 0, "value": [0.5, 0.5, 0.729167246]},
            {"id": 1, "value": [0.5, 0.0, 0.814951628]},
            {"id": 2, "value": [0.0, 0.0, 0.90073601]},
            {"id": 3, "value": [0.0, 0.5, 0.986520391]},
            {"id": 4, "value": [0.5, 0.5, 0.386029718]},
            {"id": 5, "value": [0.5, 0.0, 0.4718141]},
            {"id": 6, "value": [0.0, 0.0, 0.557598482]},
            {"id": 7, "value": [0.0, 0.5, 0.643382864]},
        ],
        "units": "crystal",
        "constraints": [],
        "labels": [],
    },
    "lattice": {
        "a": 3.867,
        "b": 3.867,
        "c": 15.937527692,
        "alpha": 90.0,
        "beta": 90.0,
        "gamma": 90.0,
        "units": {"length": "angstrom", "angle": "degree"},
        "type": "TRI",
        "vectors": {
            "a": [3.867, 0.0, 0.0],
            "b": [-0.0, 3.867, 0.0],
            "c": [0.0, 0.0, 15.937527692],
            "alat": 1,
            "units": "angstrom",
        },
    },
    "isNonPeriodic": False,
    "_id": "",
    "metadata": {
        "boundaryConditions": {"type": "pbc", "offset": 0},
        "build": {
            "termination": "Si_P4/mmm_1",
            "configuration": {
                "type": "SlabConfiguration",
                "bulk": {
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
                            {"id": 0, "value": [0.5, 0.0, 0.0]},
                            {"id": 1, "value": [0.25, 0.25, 0.75]},
                            {"id": 2, "value": [0.5, 0.5, 0.5]},
                            {"id": 3, "value": [0.25, 0.75, 0.25]},
                            {"id": 4, "value": [0.0, 0.0, 0.5]},
                            {"id": 5, "value": [0.75, 0.25, 0.25]},
                            {"id": 6, "value": [0.0, 0.5, 0.0]},
                            {"id": 7, "value": [0.75, 0.75, 0.75]},
                        ],
                        "units": "crystal",
                        "constraints": [],
                        "labels": [],
                    },
                    "lattice": {
                        "a": 5.468763846,
                        "b": 5.468763846,
                        "c": 5.468763846,
                        "alpha": 90.0,
                        "beta": 90.0,
                        "gamma": 90.0,
                        "units": {"length": "angstrom", "angle": "degree"},
                        "type": "TRI",
                        "vectors": {
                            "a": [5.468763846, 0.0, 0.0],
                            "b": [-0.0, 5.468763846, 0.0],
                            "c": [0.0, 0.0, 5.468763846],
                            "alat": 1,
                            "units": "angstrom",
                        },
                    },
                    "isNonPeriodic": False,
                    "_id": "",
                    "metadata": {"boundaryConditions": {"type": "pbc", "offset": 0}},
                    "isUpdated": True,
                },
                "miller_indices": (0, 0, 1),
                "thickness": 2,
                "vacuum": 5.0,
                "xy_supercell_matrix": [[1, 0], [0, 1]],
                "use_conventional_cell": True,
                "use_orthogonal_z": True,
                "make_primitive": True,
            },
        },
    },
    "isUpdated": True,
}

SI_SLAB: Dict[str, Any] = {
    "name": "Si8(001), termination Si_P4/mmm_1, Slab",
    "basis": {
        "elements": [{"id": 0, "value": "Si"}, {"id": 1, "value": "Si"}],
        "coordinates": [
            {"id": 0, "value": [0.583333333, 0.833333333, 0.241911889]},
            {"id": 1, "value": [0.25, 0.5, 0.145147133]},
        ],
        "units": "crystal",
        "constraints": [],
        "labels": [],
    },
    "lattice": {
        "a": 3.867,
        "b": 3.867,
        "c": 8.157392279,
        "alpha": 90.0,
        "beta": 90.0,
        "gamma": 60.0,
        "units": {"length": "angstrom", "angle": "degree"},
        "type": "TRI",
        "vectors": {
            "a": [3.867, 0.0, 0.0],
            "b": [1.9335, 3.348920236, 0.0],
            "c": [0.0, 0.0, 8.157392279],
            "alat": 1,
            "units": "angstrom",
        },
    },
    "isNonPeriodic": False,
    "_id": "",
    "metadata": {
        "boundaryConditions": {"type": "pbc", "offset": 0},
        "build": {"configuration": SI_SLAB_CONFIGURATION, "termination": "Si_P4/mmm_1"},
    },
    "isUpdated": True,
}

SI_SLAB_PASSIVATED = {
    "name": "Si8(001), termination Si_P4/mmm_1, Slab H-passivated",
    "basis": {
        "elements": [
            {"id": 0, "value": "Si"},
            {"id": 1, "value": "Si"},
            {"id": 2, "value": "H"},
            {"id": 3, "value": "H"},
        ],
        "coordinates": [
            {"id": 0, "value": [0.583333333, 0.833333333, 0.548382368]},
            {"id": 1, "value": [0.25, 0.5, 0.451617612]},
            {"id": 2, "value": [0.25, 0.5, 0.270187076]},
            {"id": 3, "value": [0.583333333, 0.833333333, 0.729812904]},
        ],
        "units": "crystal",
        "labels": [],
    },
    "lattice": {
        "a": 3.867,
        "b": 3.867,
        "c": 8.157392279,
        "alpha": 90.0,
        "beta": 90.0,
        "gamma": 60.0,
        "units": {"length": "angstrom", "angle": "degree"},
        "type": "TRI",
        "vectors": {
            "a": [3.867, 0.0, 0.0],
            "b": [1.9335, 3.348920236, 0.0],
            "c": [0.0, 0.0, 8.157392279],
            "alat": 1,
            "units": "angstrom",
        },
    },
    "isNonPeriodic": False,
    "_id": "",
    "metadata": {
        "boundaryConditions": {"type": "pbc", "offset": 0},
        "build": {
            "configuration": {
                "type": "PassivationConfiguration",
                # TODO: `basis` retains "cell" leading to a mismatch in the test
                "slab": reduce(lambda d, key: d.get(key, {}), ["basis"], SI_SLAB).pop("cell", None),
                "passivant": "H",
                "bond_length": 1.48,
                "surface": "both",
            },
            "termination": "Si_P4/mmm_1",
        },
    },
    "isUpdated": True,
}

SI_SLAB_VACUUM = copy.deepcopy(SI_SLAB)
SI_SLAB_VACUUM["basis"]["coordinates"] = [
    {"id": 0, "value": [0.583333333, 0.833333333, 0.149981861]},
    {"id": 1, "value": [0.25, 0.5, 0.089989116]},
]
# SI_SLAB_VACUUM["basis"]["cell"] = [[3.867, 0.0, 0.0], [1.9335, 3.348920236, 0.0], [0.0, 0.0, 13.157392279]]
SI_SLAB_VACUUM["lattice"]["c"] = 13.157392279
SI_SLAB_VACUUM["lattice"]["vectors"]["c"] = [0.0, 0.0, 13.157392279]

clean_material = Material.create(Material.default_config)
slab_111_config = SlabConfiguration(
    bulk=clean_material,
    miller_indices=(1, 1, 1),
    thickness=4,
    vacuum=6,
    xy_supercell_matrix=[[1, 0], [0, 1]],
    use_orthogonal_z=True,
)
t_111 = get_terminations(slab_111_config)[0]
SLAB_111 = create_slab(slab_111_config, t_111)

slab_001_config = SlabConfiguration(
    bulk=clean_material,
    miller_indices=(0, 0, 1),
    thickness=3,
    vacuum=3,
    xy_supercell_matrix=[[2, 0], [0, 1]],
    use_orthogonal_z=True,
)
t_001 = get_terminations(slab_001_config)[0]
SLAB_001 = create_slab(slab_001_config, t_001)

GRAPHENE = {
    "name": "Graphene",
    "basis": {
        "elements": [{"id": 0, "value": "C"}, {"id": 1, "value": "C"}],
        "coordinates": [{"id": 0, "value": [0, 0, 0]}, {"id": 1, "value": [0.333333, 0.666667, 0]}],
        "units": "crystal",
        "constraints": [],
    },
    "lattice": {
        "a": 2.467291,
        "b": 2.467291,
        "c": 20,
        "alpha": 90,
        "beta": 90,
        "gamma": 120,
        "units": {"length": "angstrom", "angle": "degree"},
        "type": "HEX",
        "vectors": {
            "a": [2.467291, 0, 0],
            "b": [-1.233645, 2.136737, 0],
            "c": [0, 0, 20],
            "alat": 1,
            "units": "angstrom",
        },
    },
    "isNonPeriodic": False,
}

GRAPHENE_ZIGZAG_NANORIBBON: Dict[str, Any] = {
    "name": "Graphene (Zigzag nanoribbon)",
    "basis": {
        "elements": [
            {"id": 0, "value": "C"},
            {"id": 1, "value": "C"},
            {"id": 2, "value": "C"},
            {"id": 3, "value": "C"},
            {"id": 4, "value": "C"},
            {"id": 5, "value": "C"},
            {"id": 6, "value": "C"},
            {"id": 7, "value": "C"},
            {"id": 8, "value": "C"},
            {"id": 9, "value": "C"},
            {"id": 10, "value": "C"},
            {"id": 11, "value": "C"},
            {"id": 12, "value": "C"},
            {"id": 13, "value": "C"},
            {"id": 14, "value": "C"},
            {"id": 15, "value": "C"},
        ],
        "coordinates": [
            {"id": 0, "value": [0.062499937, 0.366666681, 0.5]},
            {"id": 1, "value": [0.312499937, 0.366666681, 0.5]},
            {"id": 2, "value": [0.187500062, 0.433333286, 0.5]},
            {"id": 3, "value": [0.187499937, 0.566666695, 0.5]},
            {"id": 4, "value": [0.062500062, 0.6333333, 0.5]},
            {"id": 5, "value": [0.562499937, 0.366666681, 0.5]},
            {"id": 6, "value": [0.437500062, 0.433333286, 0.5]},
            {"id": 7, "value": [0.437499937, 0.566666695, 0.5]},
            {"id": 8, "value": [0.312500062, 0.6333333, 0.5]},
            {"id": 9, "value": [0.812499938, 0.366666681, 0.5]},
            {"id": 10, "value": [0.687500062, 0.433333286, 0.5]},
            {"id": 11, "value": [0.687499937, 0.566666695, 0.5]},
            {"id": 12, "value": [0.562500062, 0.6333333, 0.5]},
            {"id": 13, "value": [0.937500063, 0.433333286, 0.5]},
            {"id": 14, "value": [0.937499937, 0.566666695, 0.5]},
            {"id": 15, "value": [0.812500063, 0.6333333, 0.5]},
        ],
        "units": "crystal",
        "constraints": [],
        "labels": [],
    },
    "lattice": {
        "a": 9.869164,
        "b": 10.683683422,
        "c": 20.0,
        "alpha": 90.0,
        "beta": 90.0,
        "gamma": 90.0,
        "units": {"length": "angstrom", "angle": "degree"},
        "type": "TRI",
        "vectors": {
            "a": [9.869164, 0.0, 0.0],
            "b": [-0.0, 10.683683422, 0.0],
            "c": [0.0, 0.0, 20.0],
            "alat": 1,
            "units": "angstrom",
        },
    },
    "isNonPeriodic": False,
    "_id": "",
    "metadata": {
        "boundaryConditions": {"type": "pbc", "offset": 0},
        "build": {
            "configuration": {
                "material": GRAPHENE,
                "width": 2,
                "length": 4,
                "vacuum_width": 3,
                "vacuum_length": 0,
                "edge_type": "zigzag",
            }
        },
    },
    "isUpdated": True,
}

GRAPHENE_ARMCHAIR_NANORIBBON = {
    "name": "Graphene (Armchair nanoribbon)",
    "basis": {
        "elements": [
            {"id": 0, "value": "C"},
            {"id": 1, "value": "C"},
            {"id": 2, "value": "C"},
            {"id": 3, "value": "C"},
            {"id": 4, "value": "C"},
            {"id": 5, "value": "C"},
            {"id": 6, "value": "C"},
            {"id": 7, "value": "C"},
            {"id": 8, "value": "C"},
            {"id": 9, "value": "C"},
            {"id": 10, "value": "C"},
            {"id": 11, "value": "C"},
            {"id": 12, "value": "C"},
            {"id": 13, "value": "C"},
            {"id": 14, "value": "C"},
            {"id": 15, "value": "C"},
        ],
        "coordinates": [
            {"id": 0, "value": [0.958333362, 0.35000006, 0.5]},
            {"id": 1, "value": [0.791666617, 0.34999996, 0.5]},
            {"id": 2, "value": [0.958333362, 0.550000054, 0.5]},
            {"id": 3, "value": [0.791666617, 0.549999954, 0.5]},
            {"id": 4, "value": [0.708333369, 0.450000057, 0.5]},
            {"id": 5, "value": [0.541666624, 0.449999957, 0.5]},
            {"id": 6, "value": [0.458333376, 0.35000006, 0.5]},
            {"id": 7, "value": [0.291666631, 0.34999996, 0.5]},
            {"id": 8, "value": [0.708333369, 0.650000051, 0.5]},
            {"id": 9, "value": [0.541666624, 0.649999951, 0.5]},
            {"id": 10, "value": [0.458333376, 0.550000054, 0.5]},
            {"id": 11, "value": [0.291666631, 0.549999954, 0.5]},
            {"id": 12, "value": [0.208333383, 0.450000057, 0.5]},
            {"id": 13, "value": [0.041666638, 0.449999957, 0.5]},
            {"id": 14, "value": [0.208333383, 0.650000051, 0.5]},
            {"id": 15, "value": [0.041666638, 0.649999951, 0.5]},
        ],
        "units": "crystal",
        "constraints": [],
        "labels": [],
    },
    "lattice": {
        "a": 8.546946738,
        "b": 12.336455,
        "c": 20.0,
        "alpha": 90.0,
        "beta": 90.0,
        "gamma": 90.0,
        "units": {"length": "angstrom", "angle": "degree"},
        "type": "TRI",
        "vectors": {
            "a": [8.546946738, 0.0, 0.0],
            "b": [-0.0, 12.336455, 0.0],
            "c": [0.0, 0.0, 20.0],
            "alat": 1,
            "units": "angstrom",
        },
    },
    "isNonPeriodic": False,
    "_id": "",
    "metadata": {
        "boundaryConditions": {"type": "pbc", "offset": 0},
        "build": {
            "configuration": {
                "material": GRAPHENE,
                "width": 2,
                "length": 4,
                "vacuum_width": 3,
                "vacuum_length": 0,
                "edge_type": "armchair",
            }
        },
    },
    "isUpdated": True,
}

GRAPHENE_ZIGZAG_NANORIBBON_PASSIVATED = {
    "name": "Graphene (Zigzag nanoribbon) H-passivated",
    "basis": {
        "elements": [
            {"id": 0, "value": "C"},
            {"id": 1, "value": "C"},
            {"id": 2, "value": "C"},
            {"id": 3, "value": "C"},
            {"id": 4, "value": "C"},
            {"id": 5, "value": "C"},
            {"id": 6, "value": "C"},
            {"id": 7, "value": "C"},
            {"id": 8, "value": "C"},
            {"id": 9, "value": "C"},
            {"id": 10, "value": "C"},
            {"id": 11, "value": "C"},
            {"id": 12, "value": "C"},
            {"id": 13, "value": "C"},
            {"id": 14, "value": "C"},
            {"id": 15, "value": "C"},
            {"id": 16, "value": "H"},
            {"id": 17, "value": "H"},
            {"id": 18, "value": "H"},
            {"id": 19, "value": "H"},
            {"id": 20, "value": "H"},
            {"id": 21, "value": "H"},
            {"id": 22, "value": "H"},
            {"id": 23, "value": "H"},
        ],
        "coordinates": [
            {"id": 0, "value": [0.062499937, 0.366666681, 0.5]},
            {"id": 1, "value": [0.312499937, 0.366666681, 0.5]},
            {"id": 2, "value": [0.187500062, 0.433333286, 0.5]},
            {"id": 3, "value": [0.187499937, 0.566666695, 0.5]},
            {"id": 4, "value": [0.062500062, 0.6333333, 0.5]},
            {"id": 5, "value": [0.562499937, 0.366666681, 0.5]},
            {"id": 6, "value": [0.437500062, 0.433333286, 0.5]},
            {"id": 7, "value": [0.437499937, 0.566666695, 0.5]},
            {"id": 8, "value": [0.312500062, 0.6333333, 0.5]},
            {"id": 9, "value": [0.812499938, 0.366666681, 0.5]},
            {"id": 10, "value": [0.687500062, 0.433333286, 0.5]},
            {"id": 11, "value": [0.687499937, 0.566666695, 0.5]},
            {"id": 12, "value": [0.562500062, 0.6333333, 0.5]},
            {"id": 13, "value": [0.937500063, 0.433333286, 0.5]},
            {"id": 14, "value": [0.937499937, 0.566666695, 0.5]},
            {"id": 15, "value": [0.812500063, 0.6333333, 0.5]},
            {"id": 16, "value": [0.062500067, 0.228137674, 0.5]},
            {"id": 17, "value": [0.312500067, 0.228137674, 0.5]},
            {"id": 18, "value": [0.062499932, 0.771862307, 0.5]},
            {"id": 19, "value": [0.562500067, 0.228137674, 0.5]},
            {"id": 20, "value": [0.312499932, 0.771862307, 0.5]},
            {"id": 21, "value": [0.812500068, 0.228137674, 0.5]},
            {"id": 22, "value": [0.562499932, 0.771862307, 0.5]},
            {"id": 23, "value": [0.812499933, 0.771862307, 0.5]},
        ],
        "units": "crystal",
        "labels": [],
    },
    "lattice": {
        "a": 9.869164,
        "b": 10.683683422,
        "c": 20.0,
        "alpha": 90.0,
        "beta": 90.0,
        "gamma": 90.0,
        "units": {"length": "angstrom", "angle": "degree"},
        "type": "TRI",
        "vectors": {
            "a": [9.869164, 0.0, 0.0],
            "b": [-0.0, 10.683683422, 0.0],
            "c": [0.0, 0.0, 20.0],
            "alat": 1,
            "units": "angstrom",
        },
    },
    "isNonPeriodic": False,
    "_id": "",
    "metadata": {
        "boundaryConditions": {"type": "pbc", "offset": 0},
        "build": {
            "configuration": {
                "type": "PassivationConfiguration",
                # TODO: `basis` retains "cell" leading to a mismatch in the test (as above)
                "slab": reduce(lambda d, key: d.get(key, {}), ["basis"], GRAPHENE_ZIGZAG_NANORIBBON).pop("cell", None),
                "passivant": "H",
                "bond_length": 1.48,
                "surface": "both",
            }
        },
    },
    "isUpdated": True,
}

GRAPHENE_NICKEL_INTERFACE = {
    "name": "C2(001)-Ni4(111), Interface, Strain 0.105pct",
    "basis": {
        "elements": [
            {"id": 0, "value": "Ni"},
            {"id": 1, "value": "Ni"},
            {"id": 2, "value": "Ni"},
            {"id": 3, "value": "C"},
            {"id": 4, "value": "C"},
        ],
        "coordinates": [
            {"id": 0, "value": [0.666666667, 0.666666667, 0.350869517]},
            {"id": 1, "value": [1.0, 0.0, 0.425701769]},
            {"id": 2, "value": [0.333333333, 0.333333333, 0.500534021]},
            {"id": 3, "value": [0.333333333, 0.333333333, 0.611447347]},
            {"id": 4, "value": [0.666666667, 0.666666667, 0.611447347]},
        ],
        "units": "crystal",
        "constraints": [],
        "labels": [
            {"id": 0, "value": 0},
            {"id": 1, "value": 0},
            {"id": 2, "value": 0},
            {"id": 3, "value": 1},
            {"id": 4, "value": 1},
        ],
    },
    "lattice": {
        "a": 2.478974,
        "b": 2.478974,
        "c": 27.048147591,
        "alpha": 90.0,
        "beta": 90.0,
        "gamma": 60.0,
        "units": {"length": "angstrom", "angle": "degree"},
        "type": "TRI",
        "vectors": {
            "a": [2.478974, 0.0, 0.0],
            "b": [1.239487, 2.14685446, 0.0],
            "c": [0.0, 0.0, 27.048147591],
            "alat": 1,
            "units": "angstrom",
        },
    },
    "isNonPeriodic": False,
    "_id": "",
    "metadata": {
        "interface_properties": {
            "film_sl_vectors": [[2.467291, 0.0, 0.0], [1.2336455, -2.136736685, -0.0]],
            "substrate_sl_vectors": [[-1.752899326, 1.752899326, 0.0], [-1.752899326, 0.0, 1.752899326]],
            "film_vectors": [[2.467291, 0.0, 0.0], [-1.2336455, 2.136736685, 0.0]],
            "substrate_vectors": [[-1.752899326, 1.752899326, 0.0], [-1.752899326, 0.0, 1.752899326]],
            "film_transformation": [[1.0, 0.0], [0.0, 1.0]],
            "substrate_transformation": [[1.0, 0.0], [0.0, 1.0]],
            "strain": [[0.004746364, -0.0, -0.0], [-0.0, 0.004746364, -0.0], [-0.0, -0.0, 0.0]],
            "von_mises_strain": 0.003164242537164297,
            "termination": "('C_P6/mmm_2', 'Ni_R-3m_1')",
            "film_thickness": 1,
            "substrate_thickness": 3,
            "mean_abs_strain": 0.0010500000000000002,
        },
        "boundaryConditions": {"type": "pbc", "offset": 0},
        "mean_abs_strain": 0.0010500000000000002,
        "build": {
            "configuration": {
                "type": "InterfaceConfiguration",
                "film_configuration": {
                    "type": "SlabConfiguration",
                    "bulk": {**GRAPHENE, "name": "C2"},
                    "miller_indices": [0, 0, 1],
                    "thickness": 1,
                    "vacuum": 0,
                    "xy_supercell_matrix": [[1, 0], [0, 1]],
                    "use_conventional_cell": True,
                    "use_orthogonal_z": True,
                    "make_primitive": False,
                },
                "substrate_configuration": {
                    "type": "SlabConfiguration",
                    "bulk": {
                        "name": "Ni4",
                        "basis": {
                            "elements": [
                                {"id": 0, "value": "Ni"},
                                {"id": 1, "value": "Ni"},
                                {"id": 2, "value": "Ni"},
                                {"id": 3, "value": "Ni"},
                            ],
                            "coordinates": [
                                {"id": 0, "value": [0.0, 0.0, 0.0]},
                                {"id": 1, "value": [0.0, 0.5, 0.5]},
                                {"id": 2, "value": [0.5, 0.0, 0.5]},
                                {"id": 3, "value": [0.5, 0.5, 0.0]},
                            ],
                            "units": "crystal",
                            "constraints": [],
                            "labels": [],
                        },
                        "lattice": {
                            "a": 3.505798652,
                            "b": 3.505798652,
                            "c": 3.505798652,
                            "alpha": 90.0,
                            "beta": 90.0,
                            "gamma": 90.0,
                            "units": {"length": "angstrom", "angle": "degree"},
                            "type": "TRI",
                            "vectors": {
                                "a": [3.505798652, 0.0, 0.0],
                                "b": [-0.0, 3.505798652, 0.0],
                                "c": [0.0, 0.0, 3.505798652],
                                "alat": 1,
                                "units": "angstrom",
                            },
                        },
                        "isNonPeriodic": False,
                        "_id": "",
                        "metadata": {"boundaryConditions": {"type": "pbc", "offset": 0}},
                        "isUpdated": True,
                    },
                    "miller_indices": [1, 1, 1],
                    "thickness": 3,
                    "vacuum": 3,
                    "xy_supercell_matrix": [[1, 0], [0, 1]],
                    "use_conventional_cell": True,
                    "use_orthogonal_z": True,
                    "make_primitive": False,
                },
                "film_termination": "C_P6/mmm_2",
                "substrate_termination": "Ni_P6/mmm_4",
                "distance_z": 3.0,
                "vacuum": 20.0,
            }
        },
    },
    "isUpdated": True,
}
