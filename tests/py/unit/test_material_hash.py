"""
Tests that Python Material.hash / Material.calculate_hash produce values identical to
the JavaScript implementation (@mat3ra/made Material.calculateHash / scaledHash).

Expected values are cross-verified against the JS distribution:
  node -e "const {Material} = require('.../dist/js/made.js'); console.log(m.calculateHash())"
"""
import pytest
from mat3ra.made.material import Material
from unit.fixtures.bulk import BULK_Si_PRIMITIVE

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

GRAPHENE_CONFIG = {
    "name": "Graphene",
    "basis": {
        "units": "crystal",
        "elements": [{"id": 0, "value": "C"}, {"id": 1, "value": "C"}],
        "coordinates": [
            {"id": 0, "value": [0.333333, 0.666667, 0.5]},
            {"id": 1, "value": [0.666667, 0.333333, 0.5]},
        ],
        "labels": [],
        "constraints": [],
    },
    "lattice": {
        "type": "HEX",
        "a": 2.464955,
        "b": 2.464956,
        "c": 19.996729,
        "alpha": 90,
        "beta": 90,
        "gamma": 120,
        "units": {"length": "angstrom", "angle": "degree"},
    },
}

SRTIO3_CONFIG = {
    "name": "SrTiO3",
    "basis": {
        "units": "crystal",
        "elements": [
            {"id": 0, "value": "Sr"},
            {"id": 1, "value": "Ti"},
            {"id": 2, "value": "O"},
            {"id": 3, "value": "O"},
            {"id": 4, "value": "O"},
        ],
        "coordinates": [
            {"id": 0, "value": [0, 0, 0]},
            {"id": 1, "value": [0.5, 0.5, 0.5]},
            {"id": 2, "value": [0.5, 0, 0.5]},
            {"id": 3, "value": [0.5, 0.5, 0]},
            {"id": 4, "value": [0, 0.5, 0.5]},
        ],
        "labels": [],
        "constraints": [],
    },
    "lattice": {
        "type": "CUB",
        "a": 3.912701,
        "b": 3.912701,
        "c": 3.912701,
        "alpha": 90,
        "beta": 90,
        "gamma": 90,
        "units": {"length": "angstrom", "angle": "degree"},
    },
}

# ---------------------------------------------------------------------------
# Expected values – cross-verified against JS dist/js/made.js
# ---------------------------------------------------------------------------

SI_FCC_HASH = "a665723ef7429caef6ca89385fe25bae"
SI_FCC_SCALED_HASH = "87e416358789a7c435eaaba0344af51d"

GRAPHENE_HASH = "6f1e125d8985705657c0fc45601f4b99"
GRAPHENE_SCALED_HASH = "3545f3a304cea50c2dd6b685aac48c83"

SRTIO3_HASH = "a36ecf32ad03149f3103913cda94f9a7"
SRTIO3_SCALED_HASH = "d104cc0e9c857c6fdab4447a7ae13071"


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "config, expected_hash",
    [
        (BULK_Si_PRIMITIVE, SI_FCC_HASH),
        (GRAPHENE_CONFIG, GRAPHENE_HASH),
        (SRTIO3_CONFIG, SRTIO3_HASH),
    ],
    ids=["Si_FCC", "Graphene", "SrTiO3"],
)
def test_calculate_hash(config, expected_hash):
    material = Material.create(config)
    material.basis.set_labels_from_list([])
    assert material.calculate_hash() == expected_hash


@pytest.mark.parametrize(
    "config, expected_hash",
    [
        (BULK_Si_PRIMITIVE, SI_FCC_HASH),
        (GRAPHENE_CONFIG, GRAPHENE_HASH),
        (SRTIO3_CONFIG, SRTIO3_HASH),
    ],
    ids=["Si_FCC", "Graphene", "SrTiO3"],
)
def test_hash_property(config, expected_hash):
    material = Material.create(config)
    material.basis.set_labels_from_list([])
    assert material.hash == expected_hash


@pytest.mark.parametrize(
    "config, expected_scaled_hash",
    [
        (BULK_Si_PRIMITIVE, SI_FCC_SCALED_HASH),
        (GRAPHENE_CONFIG, GRAPHENE_SCALED_HASH),
        (SRTIO3_CONFIG, SRTIO3_SCALED_HASH),
    ],
    ids=["Si_FCC", "Graphene", "SrTiO3"],
)
def test_scaled_hash_property(config, expected_scaled_hash):
    material = Material.create(config)
    material.basis.set_labels_from_list([])
    assert material.scaled_hash == expected_scaled_hash

