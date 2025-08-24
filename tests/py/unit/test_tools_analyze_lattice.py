import pytest
from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.lattice import LatticeMaterialAnalyzer, get_lattice_type

from .fixtures.bulk import BULK_GRAPHITE, BULK_Hf2O_MCL, BULK_Si_CONVENTIONAL, BULK_Si_PRIMITIVE, BULK_Si_PRIMITIVIZED
from .fixtures.interface.gr_ni_111_top_hcp import GRAPHENE_NICKEL_INTERFACE_TOP_HCP
from .utils import assert_two_entities_deep_almost_equal


@pytest.mark.parametrize(
    "primitive_material_config, expected_conventional_material_config, expected_primitive_material_config",
    [(BULK_Si_PRIMITIVE, BULK_Si_CONVENTIONAL, BULK_Si_PRIMITIVIZED)],
)
def test_lattice_material_analyzer(
    primitive_material_config, expected_conventional_material_config, expected_primitive_material_config
):
    primitive_cell = Material.create(primitive_material_config)
    lattice_material_analyzer = LatticeMaterialAnalyzer(material=primitive_cell)
    conventional_cell = lattice_material_analyzer.material_with_conventional_lattice
    assert_two_entities_deep_almost_equal(conventional_cell, expected_conventional_material_config)

    primitive_cell_generated = lattice_material_analyzer.material_with_primitive_lattice
    assert_two_entities_deep_almost_equal(primitive_cell_generated, expected_primitive_material_config)


@pytest.mark.parametrize(
    "material, expected_lattice_type",
    [
        (BULK_Si_PRIMITIVE, "FCC"),
        (BULK_Si_CONVENTIONAL, "FCC"),
        (GRAPHENE_NICKEL_INTERFACE_TOP_HCP, "HEX"),
        (BULK_Hf2O_MCL, "MCL"),
        (BULK_GRAPHITE, "HEX"),
    ],
)
def test_analyze_lattice_type(material, expected_lattice_type):
    result = get_lattice_type(material)
    assert result.value == expected_lattice_type
