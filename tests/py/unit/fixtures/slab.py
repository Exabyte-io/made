import copy
from functools import reduce
from typing import Any, Dict

from .cell import SI_CONVENTIONAL_CELL

SI_SLAB_001_CONFIGURATION: Dict[str, Any] = {
    "type": "SlabConfiguration",
    "bulk": SI_CONVENTIONAL_CELL,
    "miller_indices": (0, 0, 1),
    "number_of_layers": 2,
    "vacuum": 5.0,
    "xy_supercell_matrix": [[1, 0], [0, 1]],
    "use_conventional_cell": True,
}

SI_SLAB_001_BUILD_PARAMETERS: Dict[str, Any] = {
    "min_vacuum_size": 1.0,
    "reorient_lattice": True,
    "symmetrize": True,
    "make_primitive": True,
    "use_orthogonal_c": True,
}


SI_SLAB_001_2_ATOMS: Dict[str, Any] = {
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
    },
    "isNonPeriodic": False,
    "metadata": {
        "boundaryConditions": {"type": "pbc", "offset": 0},
        "build": {"configuration": SI_SLAB_001_CONFIGURATION, "termination": "Si_P4/mmm_1"},
    },
}

SI_SLAB_001: Dict[str, Any] = {
    "name": "Si8(001), termination Si_P4/mmm_2, Slab",
    "basis": {
        "constraints": [],
        "coordinates": [
            {"id": 0, "value": [0.5, 0.5, 0.042892191]},
            {"id": 1, "value": [0.5, 0.0, 0.128676573]},
            {"id": 2, "value": [0.0, 0.0, 0.214460955]},
            {"id": 3, "value": [0.0, 0.5, 0.300245337]},
            {"id": 4, "value": [0.5, 0.5, 0.386029718]},
            {"id": 5, "value": [0.5, 0.0, 0.4718141]},
            {"id": 6, "value": [0.0, 0.0, 0.557598482]},
            {"id": 7, "value": [0.0, 0.5, 0.643382864]},
        ],
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
        "labels": [],
        "units": "crystal",
    },
    "lattice": {
        "a": 3.867,
        "alpha": 90.0,
        "b": 3.867,
        "beta": 90.0,
        "c": 15.937527692,
        "gamma": 90.0,
        "type": "TRI",
        "units": {"angle": "degree", "length": "angstrom"},
    },
    "isNonPeriodic": False,
    "metadata": {
        "boundaryConditions": {"type": "pbc", "offset": 0},
        "build": {
            "termination": "Si_P4/mmm_2",
            "configuration": {
                "type": "SlabConfiguration",
                "stack_components": [
                    {
                        "crystal": SI_CONVENTIONAL_CELL,
                        "miller_indices": [0, 0, 1],
                        "number_of_repetitions": 2,
                        "termination_top": {"chemical_elements": "Si", "space_group_symmetry_label": "P4/mmm_2"},
                        "use_conventional_cell": True,
                    },
                    {"type": "VacuumConfiguration", "direction": "z", "size": 5.0},
                ],
                "xy_supercell_matrix": [[1, 0], [0, 1]],
            },
        },
    },
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
    },
    "isNonPeriodic": False,
    "metadata": {
        "boundaryConditions": {"type": "pbc", "offset": 0},
        "build": {
            "configuration": {
                "type": "PassivationConfiguration",
                # TODO: `basis` retains "cell" leading to a mismatch in the test
                "slab": reduce(lambda d, key: d.get(key, {}), ["basis"], SI_SLAB_001_2_ATOMS).pop("cell", None),
                "passivant": "H",
                "bond_length": 1.48,
                "surface": "both",
            },
            "termination": "Si_P4/mmm_1",
        },
    },
}


SI_SLAB_001_WITH_VACUUM = copy.deepcopy(SI_SLAB_001_2_ATOMS)
SI_SLAB_001_WITH_VACUUM["lattice"]["c"] = 13.157392279
# The crystal xy coordinates are the same as SI_SLAB_001, but the z-coordinates are different
SI_SLAB_001_WITH_VACUUM["basis"]["coordinates"] = [
    {"id": 0, "value": [0.583333333, 0.833333333, 0.149981861]},
    {"id": 1, "value": [0.25, 0.5, 0.089989116]},
]

SI_SLAB_DEFAULT_PARAMETERS = {
    "basis": {
        "constraints": [],
        "coordinates": [
            {"id": 0, "value": [0.5, 0.5, 0.044191992]},
            {"id": 1, "value": [0.5, 0.0, 0.132575975]},
            {"id": 2, "value": [0.0, 0.0, 0.220959958]},
            {"id": 3, "value": [0.0, 0.5, 0.309343941]},
        ],
        "elements": [
            {"id": 0, "value": "Si"},
            {"id": 1, "value": "Si"},
            {"id": 2, "value": "Si"},
            {"id": 3, "value": "Si"},
        ],
        "labels": [],
        "units": "crystal",
    },
    "isDefault": False,
    "isNonPeriodic": False,
    "lattice": {
        "a": 3.867,
        "alpha": 90.0,
        "b": 3.867,
        "beta": 90.0,
        "c": 15.468763846,
        "gamma": 90.0,
        "type": "TRI",
        "units": {"angle": "degree", "length": "angstrom"},
    },
    "metadata": {
        "boundaryConditions": {"offset": 0, "type": "pbc"},
        "build": {
            "configuration": {
                "type": "SlabConfiguration",
                "stack_components": [
                    {
                        "crystal": SI_CONVENTIONAL_CELL,
                        "miller_indices": [0, 0, 1],
                        "number_of_repetitions": 1,
                        "termination_top": {"chemical_elements": "Si", "space_group_symmetry_label": "P4/mmm_2"},
                        "use_conventional_cell": True,
                    },
                    {"type": "VacuumConfiguration", "direction": "z", "size": 10.0},
                ],
                "xy_supercell_matrix": [[1, 0], [0, 1]],
            },
            "termination": "Si_P4/mmm_2",
        },
    },
    "name": "Si8(001), termination Si_P4/mmm_2, Slab",
}
