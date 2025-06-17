from mat3ra.made.material import Material
from mat3ra.made.tools.build.nanoribbon import NanoribbonConfiguration, create_nanoribbon
from mat3ra.made.tools.build.nanoribbon.enums import EdgeTypes
import pytest

from .fixtures.monolayer import GRAPHENE
from .fixtures.nanoribbon import GRAPHENE_ARMCHAIR_NANORIBBON, GRAPHENE_ZIGZAG_NANORIBBON
from .utils import assert_two_entities_deep_almost_equal


@pytest.mark.parametrize(
    "material_config, width, length, vacuum_width, edge_type, expected_nanoribbon",
    [
        (GRAPHENE, 2, 4, 3, EdgeTypes.zigzag, GRAPHENE_ZIGZAG_NANORIBBON),
        (GRAPHENE, 2, 4, 3, EdgeTypes.armchair, GRAPHENE_ARMCHAIR_NANORIBBON),
    ],
)
def test_build_nanoribbon(material_config, width, length, vacuum_width, edge_type, expected_nanoribbon):
    config = NanoribbonConfiguration(
        material=Material.create(material_config),
        width=width,
        length=length,
        vacuum_width=vacuum_width,
        edge_type=edge_type,
    )

    nanoribbon = create_nanoribbon(config)
    assert_two_entities_deep_almost_equal(nanoribbon, expected_nanoribbon)
