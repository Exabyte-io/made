from mat3ra.made.material import Material
from mat3ra.made.tools.build.passivation import get_unique_coordination_numbers
from mat3ra.made.tools.build.passivation.builders import (
    CoordinationBasedPassivationBuilder,
    CoordinationBasedPassivationBuilderParameters,
    SurfacePassivationBuilder,
    SurfacePassivationBuilderParameters,
)
from mat3ra.made.tools.build.passivation.configuration import PassivationConfiguration
from mat3ra.utils import assertion as assertion_utils

from .fixtures import GRAPHENE_ZIGZAG_NANORIBBON, GRAPHENE_ZIGZAG_NANORIBBON_PASSIVATED, SI_SLAB, SI_SLAB_PASSIVATED


def test_passivate_surface():
    config = PassivationConfiguration(slab=Material(SI_SLAB), passivant="H", bond_length=1.48, surface="both")
    builder = SurfacePassivationBuilder(
        build_parameters=SurfacePassivationBuilderParameters(shadowing_radius=2.5, depth=2.0)
    )
    passivated_material = builder.get_material(config)
    assertion_utils.assert_deep_almost_equal(SI_SLAB_PASSIVATED, passivated_material.to_json())


def test_get_unique_coordination_numbers():
    config = PassivationConfiguration(
        slab=Material(GRAPHENE_ZIGZAG_NANORIBBON), passivant="H", bond_length=1.48, surface="both"
    )
    unique_coordination_numbers = get_unique_coordination_numbers(config, cutoff=3.0)
    assert unique_coordination_numbers == [2, 3]


def test_passivate_coordination_based():
    config = PassivationConfiguration(slab=Material(GRAPHENE_ZIGZAG_NANORIBBON), passivant="H", bond_length=1.48)
    params = CoordinationBasedPassivationBuilderParameters(
        shadowing_radius=2.5, coordination_threshold=2, bonds_to_passivate=1
    )
    builder = CoordinationBasedPassivationBuilder(build_parameters=params)
    passivated_material = builder.get_material(config)
    assertion_utils.assert_deep_almost_equal(GRAPHENE_ZIGZAG_NANORIBBON_PASSIVATED, passivated_material.to_json())
