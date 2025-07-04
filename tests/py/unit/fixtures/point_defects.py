VACANCY_DEFECT_BULK_PRIMITIVE_Si = {
    "name": "Silicon FCC with vacancy defect",
    "lattice": {
        "a": 3.867,
        "b": 3.867,
        "c": 3.867,
        "alpha": 60.0,
        "beta": 60.0,
        "gamma": 60.0,
        "units": {"length": "angstrom", "angle": "degree"},
        "type": "FCC",
    },
    "basis": {
        "elements": [{"id": 1, "value": "Si"}],
        "coordinates": [{"id": 1, "value": [0.25, 0.25, 0.25]}],
        "units": "crystal",
        "labels": [],
        "constraints": [],
    },
}

SUBSTITUTION_DEFECT_BULK_PRIMITIVE_Si = {
    "name": "Silicon FCC with substitutional defect",
    "lattice": {
        "a": 3.867,
        "b": 3.867,
        "c": 3.867,
        "alpha": 60.0,
        "beta": 60.0,
        "gamma": 60.0,
        "units": {"length": "angstrom", "angle": "degree"},
        "type": "FCC",
    },
    "basis": {
        "elements": [{"id": 1, "value": "Si"}, {"id": 2, "value": "Ge"}],
        "coordinates": [{"id": 1, "value": [0.25, 0.25, 0.25]}, {"id": 2, "value": [0.0, 0.0, 0.0]}],
        "units": "crystal",
        "labels": [],
        "constraints": [],
    },
}

INTERSTITIAL_DEFECT_BULK_PRIMITIVE_Si = {
    "name": "Silicon FCC with interstitial defect",
    "lattice": {
        "a": 3.867,
        "b": 3.867,
        "c": 3.867,
        "alpha": 60.0,
        "beta": 60.0,
        "gamma": 60.0,
        "units": {"length": "angstrom", "angle": "degree"},
        "type": "FCC",
    },
    "basis": {
        "elements": [{"id": 0, "value": "Si"}, {"id": 1, "value": "Si"}, {"id": 2, "value": "C"}],
        "coordinates": [
            {"id": 0, "value": [0.0, 0.0, 0.0]},
            {"id": 1, "value": [0.25, 0.25, 0.25]},
            {"id": 2, "value": [0.5, 0.5, 0.5]},
        ],
        "units": "crystal",
        "labels": [],
        "constraints": [],
    },
}
