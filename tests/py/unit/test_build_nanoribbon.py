from mat3ra.made.material import Material
from mat3ra.made.tools.build.nanoribbon import NanoribbonConfiguration, build_nanoribbon
from mat3ra.utils import assertion as assertion_utils

from .fixtures import GRAPHENE, GRAPHENE_ARMCHAIR_NANORIBBON, GRAPHENE_ZIGZAG_NANORIBBON


def test_build_zigzag_nanoribbon():
    config = NanoribbonConfiguration(
        material=Material(GRAPHENE),
        width=2,
        length=4,
        vacuum_width=3,
        edge_type="zigzag",
    )

    nanoribbon = build_nanoribbon(config)
    assertion_utils.assert_deep_almost_equal(GRAPHENE_ZIGZAG_NANORIBBON, nanoribbon.to_json())


def test_build_armchair_nanoribbon():
    config = NanoribbonConfiguration(
        material=Material(GRAPHENE),
        width=2,
        length=4,
        vacuum_width=3,
        edge_type="armchair",
    )

    nanoribbon = build_nanoribbon(config)
    assertion_utils.assert_deep_almost_equal(GRAPHENE_ARMCHAIR_NANORIBBON, nanoribbon.to_json())
