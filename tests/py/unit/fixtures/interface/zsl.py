from ..monolayer import GRAPHENE

GRAPHENE_NICKEL_INTERFACE = {
    "metadata": {"build": [], "boundaryConditions": {"type": "pbc", "offset": 0}},
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
