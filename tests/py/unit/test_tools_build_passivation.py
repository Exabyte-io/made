from mat3ra.made.material import Material
from mat3ra.made.tools.build.passivation import get_unique_coordination_numbers
from mat3ra.made.tools.build.passivation.builders import (
    CoordinationBasedPassivationBuilder,
    CoordinationBasedPassivationBuilderParameters,
    SurfacePassivationBuilder,
    SurfacePassivationBuilderParameters,
)
from mat3ra.made.tools.build.passivation.configuration import PassivationConfiguration

from .fixtures.nanoribbon import GRAPHENE_ZIGZAG_NANORIBBON, GRAPHENE_ZIGZAG_NANORIBBON_PASSIVATED
from .fixtures.slab import SI_SLAB_001, SI_SLAB_PASSIVATED
from .utils import assert_two_entities_deep_almost_equal


def test_passivate_surface():
    config = PassivationConfiguration(
        slab=Material.create(SI_SLAB_001), passivant="H", bond_length=1.48, surface="both"
    )
    builder = SurfacePassivationBuilder(
        build_parameters=SurfacePassivationBuilderParameters(shadowing_radius=2.5, depth=2.0)
    )
    passivated_material = builder.get_material(config)
    assert_two_entities_deep_almost_equal(passivated_material, SI_SLAB_PASSIVATED)


def test_get_unique_coordination_numbers():
    config = PassivationConfiguration(
        slab=Material.create(GRAPHENE_ZIGZAG_NANORIBBON), passivant="H", bond_length=1.48, surface="both"
    )
    unique_coordination_numbers = get_unique_coordination_numbers(config, cutoff=3.0)
    assert unique_coordination_numbers == [2, 3]


def test_passivate_coordination_based():
    config = PassivationConfiguration(slab=Material.create(GRAPHENE_ZIGZAG_NANORIBBON), passivant="H", bond_length=1.48)
    params = CoordinationBasedPassivationBuilderParameters(
        shadowing_radius=2.5, coordination_threshold=2, bonds_to_passivate=1
    )
    builder = CoordinationBasedPassivationBuilder(build_parameters=params)
    passivated_material = builder.get_material(config)
    assert_two_entities_deep_almost_equal(passivated_material, GRAPHENE_ZIGZAG_NANORIBBON_PASSIVATED)
