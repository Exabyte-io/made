import pytest

from mat3ra.made.material import Material
from mat3ra.made.tools.build import MaterialWithBuildMetadata
from mat3ra.made.tools.build.passivation.configuration import PassivationConfiguration
from mat3ra.made.tools.build.passivation.helpers import (
    create_passivated_surface,
    passivate_dangling_bonds,
    get_unique_coordination_numbers,
    get_coordination_numbers_distribution,
)
from .fixtures.nanoribbon.nanoribbon import GRAPHENE_ZIGZAG_NANORIBBON, GRAPHENE_ZIGZAG_NANORIBBON_PASSIVATED
from .fixtures.slab import SI_SLAB_PASSIVATED, SI_CONVENTIONAL_SLAB_001
from .utils import assert_two_entities_deep_almost_equal


@pytest.mark.parametrize(
    "slab_config, passivant, bond_length, surface, passivation_parameters, expected_material_config",
    [
        (
            SI_CONVENTIONAL_SLAB_001,
            "H",
            1.48,
            "both",
            {"shadowing_radius": 2.5, "depth": 2.0},
            SI_SLAB_PASSIVATED,
        ),
    ],
)
def test_passivate_surface(
    slab_config, passivant, bond_length, surface, passivation_parameters, expected_material_config
):
    slab = MaterialWithBuildMetadata.create(slab_config)
    passivated_material = create_passivated_surface(slab, passivant, bond_length, **passivation_parameters)
    assert_two_entities_deep_almost_equal(passivated_material, expected_material_config)


#
@pytest.mark.parametrize(
    "slab_config, passivant, bond_length, surface, cutoff, expected_numbers",
    [(GRAPHENE_ZIGZAG_NANORIBBON, "H", 1.48, "both", 3.0, [2, 3])],
)
def test_get_unique_coordination_numbers(slab_config, passivant, bond_length, surface, cutoff, expected_numbers):
    config = PassivationConfiguration(
        material=Material.create(slab_config), passivant=passivant, bond_length=bond_length, surface=surface
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
        material=Material.create(slab_config), passivant=passivant, bond_length=bond_length, surface=surface
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
            {"shadowing_radius": 2.5, "coordination_threshold": 2, "number_of_bonds_to_passivate": 1},
            GRAPHENE_ZIGZAG_NANORIBBON_PASSIVATED,
        ),
    ],
)
def test_passivate_coordination_based(
    slab_config, passivant, bond_length, builder_params_dict, expected_material_config
):
    passivated_material = passivate_dangling_bonds(
        MaterialWithBuildMetadata.create(slab_config),
        passivant=passivant,
        bond_length=bond_length,
        **builder_params_dict,
    )
    assert_two_entities_deep_almost_equal(passivated_material, expected_material_config)
