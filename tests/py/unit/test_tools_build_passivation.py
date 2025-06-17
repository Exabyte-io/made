import pytest
from mat3ra.made.material import Material
from mat3ra.made.tools.build.passivation import get_coordination_numbers_distribution, get_unique_coordination_numbers
from mat3ra.made.tools.build.passivation.builders import (
    CoordinationBasedPassivationBuilder,
    CoordinationBasedPassivationBuilderParameters,
    SurfacePassivationBuilder,
    SurfacePassivationBuilderParameters,
)
from mat3ra.made.tools.build.passivation.configuration import PassivationConfiguration

from .fixtures.nanoribbon import GRAPHENE_ZIGZAG_NANORIBBON, GRAPHENE_ZIGZAG_NANORIBBON_PASSIVATED
from .fixtures.slab import SI_SLAB_001_2_ATOMS, SI_SLAB_PASSIVATED
from .utils import assert_two_entities_deep_almost_equal


@pytest.mark.parametrize(
    "slab_config, passivant, bond_length, surface, builder_params_dict, expected_material_config",
    [
        (
            SI_SLAB_001_2_ATOMS,
            "H",
            1.48,
            "both",
            {"shadowing_radius": 2.5, "depth": 2.0},
            SI_SLAB_PASSIVATED,
        ),
    ],
)
def test_passivate_surface(slab_config, passivant, bond_length, surface, builder_params_dict, expected_material_config):
    # TODO: use Silicon SLAB with vacuum same as for adatom
    config = PassivationConfiguration(
        slab=Material.create(slab_config), passivant=passivant, bond_length=bond_length, surface=surface
    )
    builder = SurfacePassivationBuilder(build_parameters=SurfacePassivationBuilderParameters(**builder_params_dict))
    passivated_material = builder.get_material(config)
    assert_two_entities_deep_almost_equal(passivated_material, expected_material_config)


@pytest.mark.parametrize(
    "slab_config, passivant, bond_length, surface, cutoff, expected_numbers",
    [(GRAPHENE_ZIGZAG_NANORIBBON, "H", 1.48, "both", 3.0, [2, 3])],
)
def test_get_unique_coordination_numbers(slab_config, passivant, bond_length, surface, cutoff, expected_numbers):
    config = PassivationConfiguration(
        slab=Material.create(slab_config), passivant=passivant, bond_length=bond_length, surface=surface
    )
    unique_coordination_numbers = get_unique_coordination_numbers(config, cutoff=cutoff)
    assert unique_coordination_numbers == expected_numbers


@pytest.mark.parametrize(
    "slab_config, passivant, bond_length, surface, cutoff, expected_distribution",
    [(GRAPHENE_ZIGZAG_NANORIBBON, "H", 1.48, "both", 3.0, {2: 8, 3: 8})],
)
def test_get_coordination_numbers_distribution(
    slab_config, passivant, bond_length, surface, cutoff, expected_distribution
):
    """Test getting coordination numbers distribution for passivation analysis"""
    config = PassivationConfiguration(
        slab=Material.create(slab_config), passivant=passivant, bond_length=bond_length, surface=surface
    )
    distribution = get_coordination_numbers_distribution(config, cutoff=cutoff)
    # Should return a dictionary with coordination numbers as keys and counts as values
    assert distribution == expected_distribution


@pytest.mark.parametrize(
    "slab_config, passivant, bond_length, builder_params_dict, expected_material_config",
    [
        (
            GRAPHENE_ZIGZAG_NANORIBBON,
            "H",
            1.48,
            {"shadowing_radius": 2.5, "coordination_threshold": 2, "bonds_to_passivate": 1},
            GRAPHENE_ZIGZAG_NANORIBBON_PASSIVATED,
        ),
    ],
)
def test_passivate_coordination_based(
    slab_config, passivant, bond_length, builder_params_dict, expected_material_config
):
    config = PassivationConfiguration(slab=Material.create(slab_config), passivant=passivant, bond_length=bond_length)
    params = CoordinationBasedPassivationBuilderParameters(**builder_params_dict)
    builder = CoordinationBasedPassivationBuilder(build_parameters=params)
    passivated_material = builder.get_material(config)
    assert_two_entities_deep_almost_equal(passivated_material, expected_material_config)
