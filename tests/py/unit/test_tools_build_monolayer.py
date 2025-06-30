import pytest
from mat3ra.made.material import Material
from mat3ra.made.tools.build.monolayer.helpers import create_monolayer

from .fixtures.bulk import BULK_Si_PRIMITIVE
from .fixtures.generated.fixtures import BULK_GRAPHITE
from .fixtures.monolayer import GRAPHENE, SILICENE
from .utils import assert_two_entities_deep_almost_equal


@pytest.mark.parametrize(
    "material_config, vacuum, expected_material_config",
    [
        (BULK_Si_PRIMITIVE, 10.0, SILICENE),
        (BULK_GRAPHITE, 12.197, GRAPHENE),
    ],
)
def test_create_monolayer(material_config, vacuum, expected_material_config):
    crystal = Material.create(material_config)
    monolayer = create_monolayer(crystal, vacuum=vacuum)

    assert isinstance(monolayer, Material)
    assert "Monolayer" in monolayer.name

    print("Monolayer material:", monolayer.to_dict())

    # reset the name and lattice.type to ignore them in the comparison
    monolayer.name = expected_material_config["name"]
    monolayer.lattice.type = expected_material_config["lattice"]["type"]
    monolayer.metadata.pop("build")

    assert_two_entities_deep_almost_equal(monolayer, expected_material_config, atol=1e-5)
