import pytest

from mat3ra.made.material import Material
from mat3ra.made.tools.build.nanoparticle.helpers import create_nanoparticle_by_condition
from mat3ra.made.tools.utils.coordinate import SphereCoordinateCondition
from unit.fixtures.bulk import BULK_Si_PRIMITIVE


@pytest.mark.parametrize(
    "material_config, condition, expected_material_config",
    [
        (
            BULK_Si_PRIMITIVE,
            SphereCoordinateCondition(radius=6.0),
            BULK_Si_PRIMITIVE,  # Change to nanoparticle material
        ),
    ],
)
def test_create_nanoparticle(material_config, condition, expected_material_config):
    material = Material.create(material_config)
    nanoparticle = create_nanoparticle_by_condition(
        material, condition=condition, center_around_atom=True, use_cartesian_coordinates=True
    )
    assert isinstance(nanoparticle, Material)
