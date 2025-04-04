from mat3ra.made.material import Material
from mat3ra.made.tools.build.nanoribbon import NanoribbonConfiguration, create_nanoribbon
from mat3ra.made.tools.build.nanoribbon.enums import EdgeTypes

from .fixtures.monolayer import GRAPHENE
from .fixtures.nanoribbon import GRAPHENE_ARMCHAIR_NANORIBBON, GRAPHENE_ZIGZAG_NANORIBBON
from .utils import assert_two_entities_deep_almost_equal


def test_build_zigzag_nanoribbon():
    config = NanoribbonConfiguration(
        material=Material.create(GRAPHENE),
        width=2,
        length=4,
        vacuum_width=3,
        edge_type=EdgeTypes.zigzag,
    )

    nanoribbon = create_nanoribbon(config)
    assert_two_entities_deep_almost_equal(nanoribbon, GRAPHENE_ZIGZAG_NANORIBBON)


def test_build_armchair_nanoribbon():
    config = NanoribbonConfiguration(
        material=Material.create(GRAPHENE),
        width=2,
        length=4,
        vacuum_width=3,
        edge_type=EdgeTypes.armchair,
    )

    nanoribbon = create_nanoribbon(config)
    assert_two_entities_deep_almost_equal(nanoribbon, GRAPHENE_ARMCHAIR_NANORIBBON)
