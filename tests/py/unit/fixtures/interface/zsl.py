# ruff: noqa
GRAPHENE_NICKEL_INTERFACE = {
    "metadata": {
        "build": [
            {
                "configuration": {
                    "stack_components": [
                        {
                            "xy_supercell_matrix": [[1, 0], [0, 1]],
                            "strain_matrix": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                            "stack_components": [
                                {
                                    "miller_indices": [1, 1, 1],
                                    "crystal": {
                                        "metadata": {"build": []},
                                        "name": "Ni, Nickel, FCC (Fm-3m) 3D (Bulk), mp-23",
                                        "isDefault": False,
                                        "basis": {
                                            "elements": [{"id": 0, "value": "Ni"}],
                                            "coordinates": [{"id": 0, "value": [0, 0, 0]}],
                                            "units": "crystal",
                                            "labels": [],
                                            "constraints": [],
                                        },
                                        "lattice": {
                                            "a": 2.478974,
                                            "b": 2.478974,
                                            "c": 2.478974,
                                            "alpha": 60,
                                            "beta": 60,
                                            "gamma": 60,
                                            "units": {"length": "angstrom", "angle": "degree"},
                                            "type": "FCC",
                                        },
                                        "external": {
                                            "id": "mp-23",
                                            "source": "Materials Project",
                                            "origin": True,
                                            "doi": "10.17188/1199153",
                                            "url": "https://next-gen.materialsproject.org/materials/mp-23",
                                        },
                                        "isNonPeriodic": False,
                                    },
                                    "use_conventional_cell": True,
                                    "termination_top": {
                                        "chemical_elements": "Ni",
                                        "space_group_symmetry_label": "R-3m_4",
                                    },
                                    "number_of_repetitions": 3,
                                },
                                {
                                    "direction": "z",
                                    "size": 0,
                                    "crystal": {
                                        "metadata": {
                                            "build": [
                                                {
                                                    "configuration": {
                                                        "miller_indices": [1, 1, 1],
                                                        "crystal": {
                                                            "metadata": {"build": []},
                                                            "name": "Ni, Nickel, FCC (Fm-3m) 3D (Bulk), mp-23",
                                                            "isDefault": False,
                                                            "basis": {
                                                                "elements": [{"id": 0, "value": "Ni"}],
                                                                "coordinates": [{"id": 0, "value": [0, 0, 0]}],
                                                                "units": "crystal",
                                                                "labels": [],
                                                                "constraints": [],
                                                            },
                                                            "lattice": {
                                                                "a": 2.478974,
                                                                "b": 2.478974,
                                                                "c": 2.478974,
                                                                "alpha": 60,
                                                                "beta": 60,
                                                                "gamma": 60,
                                                                "units": {"length": "angstrom", "angle": "degree"},
                                                                "type": "FCC",
                                                            },
                                                            "external": {
                                                                "id": "mp-23",
                                                                "source": "Materials Project",
                                                                "origin": True,
                                                                "doi": "10.17188/1199153",
                                                                "url": "https://next-gen.materialsproject.org/materials/mp-23",
                                                            },
                                                            "isNonPeriodic": False,
                                                        },
                                                        "use_conventional_cell": True,
                                                        "termination_top": {
                                                            "chemical_elements": "Ni",
                                                            "space_group_symmetry_label": "R-3m_4",
                                                        },
                                                        "number_of_repetitions": 3,
                                                    },
                                                    "build_parameters": {},
                                                }
                                            ],
                                            "boundaryConditions": {"type": "pbc", "offset": 0},
                                        },
                                        "name": "Ni(111), termination Ni_R-3m_4",
                                        "isDefault": False,
                                        "formula": "Ni",
                                        "basis": {
                                            "elements": [
                                                {"id": 0, "value": "Ni"},
                                                {"id": 1, "value": "Ni"},
                                                {"id": 2, "value": "Ni"},
                                            ],
                                            "coordinates": [
                                                {"id": 0, "value": [0, 0, 0.000001]},
                                                {"id": 1, "value": [0, 0, 0.333334333]},
                                                {"id": 2, "value": [0, 0, 0.666667667]},
                                            ],
                                            "units": "crystal",
                                            "labels": [],
                                            "constraints": [],
                                        },
                                        "lattice": {
                                            "a": 2.478974,
                                            "b": 2.478974,
                                            "c": 7.436922,
                                            "alpha": 120,
                                            "beta": 120,
                                            "gamma": 60,
                                            "units": {"length": "angstrom", "angle": "degree"},
                                            "type": "TRI",
                                        },
                                        "isNonPeriodic": False,
                                    },
                                    "type": "VacuumConfiguration",
                                },
                            ],
                            "direction": "z",
                            "type": "SlabStrainedSupercellConfiguration",
                            "gaps": [],
                        },
                        {
                            "xy_supercell_matrix": [[1, 0], [0, 1]],
                            "strain_matrix": [
                                [1.004735152845773, 0, 0],
                                [1.1601682219195737, 1.0047351528457729, 0],
                                [0, 0, 1],
                            ],
                            "stack_components": [
                                {
                                    "miller_indices": [0, 0, 1],
                                    "crystal": {
                                        "metadata": {"build": []},
                                        "name": "Graphene",
                                        "isDefault": False,
                                        "basis": {
                                            "elements": [{"id": 0, "value": "C"}, {"id": 1, "value": "C"}],
                                            "coordinates": [
                                                {"id": 0, "value": [0, 0, 0]},
                                                {"id": 1, "value": [0.333333, 0.666667, 0]},
                                            ],
                                            "units": "crystal",
                                            "labels": [],
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
                                        },
                                        "isNonPeriodic": False,
                                    },
                                    "use_conventional_cell": True,
                                    "termination_top": {
                                        "chemical_elements": "C",
                                        "space_group_symmetry_label": "P6/mmm_2",
                                    },
                                    "number_of_repetitions": 1,
                                },
                                {
                                    "direction": "z",
                                    "size": 0,
                                    "crystal": {
                                        "metadata": {
                                            "build": [
                                                {
                                                    "configuration": {
                                                        "miller_indices": [0, 0, 1],
                                                        "crystal": {
                                                            "metadata": {"build": []},
                                                            "name": "Graphene",
                                                            "isDefault": False,
                                                            "basis": {
                                                                "elements": [
                                                                    {"id": 0, "value": "C"},
                                                                    {"id": 1, "value": "C"},
                                                                ],
                                                                "coordinates": [
                                                                    {"id": 0, "value": [0, 0, 0]},
                                                                    {"id": 1, "value": [0.333333, 0.666667, 0]},
                                                                ],
                                                                "units": "crystal",
                                                                "labels": [],
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
                                                            },
                                                            "isNonPeriodic": False,
                                                        },
                                                        "use_conventional_cell": True,
                                                        "termination_top": {
                                                            "chemical_elements": "C",
                                                            "space_group_symmetry_label": "P6/mmm_2",
                                                        },
                                                        "number_of_repetitions": 1,
                                                    },
                                                    "build_parameters": {},
                                                }
                                            ],
                                            "boundaryConditions": {"type": "pbc", "offset": 0},
                                        },
                                        "name": "C(001), termination C_P6/mmm_2",
                                        "isDefault": False,
                                        "formula": "C",
                                        "basis": {
                                            "elements": [{"id": 0, "value": "C"}, {"id": 1, "value": "C"}],
                                            "coordinates": [
                                                {"id": 0, "value": [0, 0, 0.000001]},
                                                {"id": 1, "value": [0.333333, 0.666667, 0.000001]},
                                            ],
                                            "units": "crystal",
                                            "labels": [],
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
                                            "type": "TRI",
                                        },
                                        "isNonPeriodic": False,
                                    },
                                    "type": "VacuumConfiguration",
                                },
                            ],
                            "direction": "z",
                            "type": "SlabStrainedSupercellConfiguration",
                            "gaps": [],
                        },
                        {
                            "direction": "z",
                            "size": 10,
                            "crystal": {
                                "metadata": {
                                    "build": [
                                        {
                                            "configuration": {
                                                "miller_indices": [1, 1, 1],
                                                "crystal": {
                                                    "metadata": {"build": []},
                                                    "name": "Ni, Nickel, FCC (Fm-3m) 3D (Bulk), mp-23",
                                                    "isDefault": False,
                                                    "basis": {
                                                        "elements": [{"id": 0, "value": "Ni"}],
                                                        "coordinates": [{"id": 0, "value": [0, 0, 0]}],
                                                        "units": "crystal",
                                                        "labels": [],
                                                        "constraints": [],
                                                    },
                                                    "lattice": {
                                                        "a": 2.478974,
                                                        "b": 2.478974,
                                                        "c": 2.478974,
                                                        "alpha": 60,
                                                        "beta": 60,
                                                        "gamma": 60,
                                                        "units": {"length": "angstrom", "angle": "degree"},
                                                        "type": "FCC",
                                                    },
                                                    "external": {
                                                        "id": "mp-23",
                                                        "source": "Materials Project",
                                                        "origin": True,
                                                        "doi": "10.17188/1199153",
                                                        "url": "https://next-gen.materialsproject.org/materials/mp-23",
                                                    },
                                                    "isNonPeriodic": False,
                                                },
                                                "use_conventional_cell": True,
                                                "termination_top": {
                                                    "chemical_elements": "Ni",
                                                    "space_group_symmetry_label": "R-3m_4",
                                                },
                                                "number_of_repetitions": 3,
                                            },
                                            "build_parameters": {},
                                        },
                                        {
                                            "configuration": {
                                                "xy_supercell_matrix": [[1, 0], [0, 1]],
                                                "strain_matrix": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                                                "stack_components": [
                                                    {
                                                        "miller_indices": [1, 1, 1],
                                                        "crystal": {
                                                            "metadata": {"build": []},
                                                            "name": "Ni, Nickel, FCC (Fm-3m) 3D (Bulk), mp-23",
                                                            "isDefault": False,
                                                            "basis": {
                                                                "elements": [{"id": 0, "value": "Ni"}],
                                                                "coordinates": [{"id": 0, "value": [0, 0, 0]}],
                                                                "units": "crystal",
                                                                "labels": [],
                                                                "constraints": [],
                                                            },
                                                            "lattice": {
                                                                "a": 2.478974,
                                                                "b": 2.478974,
                                                                "c": 2.478974,
                                                                "alpha": 60,
                                                                "beta": 60,
                                                                "gamma": 60,
                                                                "units": {"length": "angstrom", "angle": "degree"},
                                                                "type": "FCC",
                                                            },
                                                            "external": {
                                                                "id": "mp-23",
                                                                "source": "Materials Project",
                                                                "origin": True,
                                                                "doi": "10.17188/1199153",
                                                                "url": "https://next-gen.materialsproject.org/materials/mp-23",
                                                            },
                                                            "isNonPeriodic": False,
                                                        },
                                                        "use_conventional_cell": True,
                                                        "termination_top": {
                                                            "chemical_elements": "Ni",
                                                            "space_group_symmetry_label": "R-3m_4",
                                                        },
                                                        "number_of_repetitions": 3,
                                                    },
                                                    {
                                                        "direction": "z",
                                                        "size": 0,
                                                        "crystal": {
                                                            "metadata": {
                                                                "build": [
                                                                    {
                                                                        "configuration": {
                                                                            "miller_indices": [1, 1, 1],
                                                                            "crystal": {
                                                                                "metadata": {"build": []},
                                                                                "name": "Ni, Nickel, FCC (Fm-3m) 3D (Bulk), mp-23",
                                                                                "isDefault": False,
                                                                                "basis": {
                                                                                    "elements": [
                                                                                        {"id": 0, "value": "Ni"}
                                                                                    ],
                                                                                    "coordinates": [
                                                                                        {"id": 0, "value": [0, 0, 0]}
                                                                                    ],
                                                                                    "units": "crystal",
                                                                                    "labels": [],
                                                                                    "constraints": [],
                                                                                },
                                                                                "lattice": {
                                                                                    "a": 2.478974,
                                                                                    "b": 2.478974,
                                                                                    "c": 2.478974,
                                                                                    "alpha": 60,
                                                                                    "beta": 60,
                                                                                    "gamma": 60,
                                                                                    "units": {
                                                                                        "length": "angstrom",
                                                                                        "angle": "degree",
                                                                                    },
                                                                                    "type": "FCC",
                                                                                },
                                                                                "derivedProperties": None,
                                                                                "external": {
                                                                                    "id": "mp-23",
                                                                                    "source": "Materials Project",
                                                                                    "origin": True,
                                                                                    "data": None,
                                                                                    "doi": "10.17188/1199153",
                                                                                    "url": "https://next-gen.materialsproject.org/materials/mp-23",
                                                                                },
                                                                            },
                                                                            "use_conventional_cell": True,
                                                                            "termination_top": {
                                                                                "chemical_elements": "Ni",
                                                                                "space_group_symmetry_label": "R-3m_4",
                                                                            },
                                                                            "termination_bottom": None,
                                                                            "number_of_repetitions": 3,
                                                                        },
                                                                        "build_parameters": {},
                                                                    }
                                                                ],
                                                                "boundaryConditions": {"type": "pbc", "offset": 0},
                                                            },
                                                            "name": "Ni(111), termination Ni_R-3m_4",
                                                            "formula": "Ni",
                                                            "basis": {
                                                                "elements": [
                                                                    {"id": 0, "value": "Ni"},
                                                                    {"id": 1, "value": "Ni"},
                                                                    {"id": 2, "value": "Ni"},
                                                                ],
                                                                "coordinates": [
                                                                    {"id": 0, "value": [0, 0, 0.000001]},
                                                                    {"id": 1, "value": [0, 0, 0.333334333]},
                                                                    {"id": 2, "value": [0, 0, 0.666667667]},
                                                                ],
                                                                "units": "crystal",
                                                                "labels": [],
                                                                "constraints": [],
                                                            },
                                                            "lattice": {
                                                                "a": 2.478974,
                                                                "b": 2.478974,
                                                                "c": 7.436922,
                                                                "alpha": 120,
                                                                "beta": 120,
                                                                "gamma": 60,
                                                                "units": {"length": "angstrom", "angle": "degree"},
                                                                "type": "TRI",
                                                            },
                                                        },
                                                        "type": "VacuumConfiguration",
                                                    },
                                                ],
                                                "direction": "z",
                                                "type": "SlabStrainedSupercellConfiguration",
                                                "gaps": [],
                                            },
                                            "build_parameters": {
                                                "use_orthogonal_c": True,
                                                "xy_supercell_matrix": [[1, 0], [0, 1]],
                                            },
                                        },
                                    ],
                                    "boundaryConditions": {"type": "pbc", "offset": 0},
                                },
                                "name": "Ni(111), termination Ni_R-3m_4, Slab",
                                "formula": "Ni",
                                "basis": {
                                    "elements": [
                                        {"id": 0, "value": "Ni"},
                                        {"id": 1, "value": "Ni"},
                                        {"id": 2, "value": "Ni"},
                                    ],
                                    "coordinates": [
                                        {"id": 0, "value": [-0.000001, -0.000001, 0.000001]},
                                        {"id": 1, "value": [0.666665667, 0.666665667, 0.333334333]},
                                        {"id": 2, "value": [0.333332333, 0.333332333, 0.666667667]},
                                    ],
                                    "units": "crystal",
                                    "labels": [
                                        {"id": 0, "value": "0"},
                                        {"id": 1, "value": "0"},
                                        {"id": 2, "value": "0"},
                                    ],
                                    "constraints": [],
                                },
                                "lattice": {
                                    "a": 2.478974,
                                    "b": 2.478974,
                                    "c": 6.072221386,
                                    "alpha": 90,
                                    "beta": 90,
                                    "gamma": 59.999999999999986,
                                    "units": {"length": "angstrom", "angle": "degree"},
                                    "type": "TRI",
                                },
                            },
                            "type": "VacuumConfiguration",
                        },
                    ],
                    "direction": "z",
                    "xy_shift": [0, 0],
                    "type": "InterfaceConfiguration",
                    "gaps": [{"id": 0, "value": 3}, {"id": 1, "value": 3}],
                },
                "build_parameters": {},
            },
        ],
        "boundaryConditions": {"type": "pbc", "offset": 0},
    },
    "name": "C(001)-Ni(111), Interface, Strain 0.335pct",
    "basis": {
        "elements": [
            {"id": 0, "value": "Ni"},
            {"id": 1, "value": "Ni"},
            {"id": 2, "value": "Ni"},
            {"id": 3, "value": "C"},
            {"id": 4, "value": "C"},
        ],
        "coordinates": [
            {"id": 0, "value": [0.999999, 0.999999, 3.03e-7]},
            {"id": 1, "value": [0.666665667, 0.666665667, 0.100960811]},
            {"id": 2, "value": [0.333332333, 0.333332333, 0.201921319]},
            {"id": 3, "value": [0, 0, 0.351561882]},
            {"id": 4, "value": [0.666666, 0.666667, 0.351561882]},
        ],
        "units": "crystal",
        "labels": [
            {"id": 0, "value": 0},
            {"id": 1, "value": 0},
            {"id": 2, "value": 0},
            {"id": 3, "value": 1},
            {"id": 4, "value": 1},
        ],
        "constraints": [],
    },
    "lattice": {
        "a": 2.478974,
        "b": 2.478974,
        "c": 20.048173665,
        "alpha": 90,
        "beta": 90,
        "gamma": 60,
        "units": {"length": "angstrom", "angle": "degree"},
        "type": "TRI",
    },
}
