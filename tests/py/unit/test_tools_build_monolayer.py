import pytest

from mat3ra.made.material import Material
from mat3ra.made.tools.build.monolayer.helpers import create_monolayer
from .fixtures.bulk import BULK_Si_PRIMITIVE
from .fixtures.generated.fixtures import BULK_GRAPHITE, SILICENE
from .fixtures.monolayer import GRAPHENE
from .utils import assert_two_entities_deep_almost_equal, assert_slab_structures_almost_equal

BULK_Si_PRIMITIVE_TEMP = {
    "name": "Silicon FCC",
    "basis": {
        "constraints": [],
        "coordinates": [{"id": 0, "value": [0.0, 0.0, 0.0]}, {"id": 1, "value": [0.55, 0.25, 0.25]}],
        "elements": [{"id": 0, "value": "Si"}, {"id": 1, "value": "Si"}],
        "labels": [],
        "units": "crystal",
    },
    "lattice": {
        "a": 3.867,
        "alpha": 60.0,
        "b": 3.867,
        "beta": 60.0,
        "c": 3.867,
        "gamma": 60.0,
        "type": "FCC",
        "units": {"angle": "degree", "length": "angstrom"},
    },
}


@pytest.mark.parametrize(
    "material_config, vacuum, expected_material_config",
    [
        (BULK_Si_PRIMITIVE, 10.0, SILICENE),
        (BULK_Si_PRIMITIVE_TEMP, 10.0, SILICENE),
        (BULK_GRAPHITE, 12.197, GRAPHENE),
    ],
)
def test_create_monolayer(material_config, vacuum, expected_material_config):
    crystal = Material.create(material_config)
    monolayer = create_monolayer(crystal, vacuum=vacuum)

    assert isinstance(monolayer, Material)
    assert "Monolayer" in monolayer.name
    monolayer.metadata.pop("build")
    monolayer.name = expected_material_config["name"]
    monolayer.lattice.type = expected_material_config["lattice"]["type"]

    assert_slab_structures_almost_equal(monolayer, expected_material_config)
