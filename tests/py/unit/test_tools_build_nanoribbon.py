import pytest
from mat3ra.made.material import Material
from mat3ra.made.tools.build.nanoribbon import (
    NanoribbonConfiguration,
    create_nanoribbon,
    CrystalLatticeLinesConfiguration,
)
from mat3ra.made.tools.build.nanoribbon.builders import NanoribbonBuilder, CrystalLatticeLinesBuilder

from .fixtures.monolayer import GRAPHENE
from .fixtures.nanoribbon import GRAPHENE_ARMCHAIR_NANORIBBON, GRAPHENE_ZIGZAG_NANORIBBON
from .utils import assert_two_entities_deep_almost_equal


@pytest.mark.parametrize(
    "material_config, width, length, vacuum_width, vacuum_length, miller_indices_uv, edge_type_name, expected_nanoribbon",
    [
        (GRAPHENE, 2, 4, 10.0, 0.0, (1, 1), "zigzag", GRAPHENE_ZIGZAG_NANORIBBON),
        (GRAPHENE, 2, 4, 10.0, 0.0, (0, 1), "armchair", GRAPHENE_ARMCHAIR_NANORIBBON),
    ],
)
def test_build_nanoribbon(
    material_config, width, length, vacuum_width, vacuum_length, miller_indices_uv, edge_type_name, expected_nanoribbon
):
    """Test nanoribbon creation using Miller indices."""
    material = Material.create(material_config)

    # Test the main API function
    nanoribbon = create_nanoribbon(
        material=material,
        miller_indices_uv=miller_indices_uv,
        width=width,
        length=length,
        vacuum_width=vacuum_width,
        vacuum_length=vacuum_length,
    )

    assert_two_entities_deep_almost_equal(nanoribbon, expected_nanoribbon)


def test_crystal_lattice_lines():
    builder = CrystalLatticeLinesBuilder()
    config = CrystalLatticeLinesConfiguration(crystal=Material.create(GRAPHENE), miller_indices_uv=(1, 0))

    material = builder.get_material(config)
    material.name = f"{material.name} {config.miller_indices_uv}"
    assert isinstance(material, Material)
