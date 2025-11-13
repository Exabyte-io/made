import pytest

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.lattice_swap_analyzer import MaterialLatticeSwapAnalyzer
from .fixtures.slab import (
    SLAB_SrTiO3_011_TERMINATION_O2,
    SLAB_SrTiO3_011_TERMINATION_SrTiO,
)
from .utils import assert_two_entities_deep_almost_equal


@pytest.mark.parametrize(
    "original_material_config, another_material_config, is_swapped",
    [
        (SLAB_SrTiO3_011_TERMINATION_O2, SLAB_SrTiO3_011_TERMINATION_O2, False),
        (SLAB_SrTiO3_011_TERMINATION_O2, SLAB_SrTiO3_011_TERMINATION_SrTiO, True),
    ],
)
def test_flip_detection_between_materials(original_material_config, another_material_config, is_swapped):
    """Test rotation detection between different materials."""
    original_material = Material.create(original_material_config)
    another_material = Material.create(another_material_config)

    rotation_analyzer = MaterialLatticeSwapAnalyzer(material=another_material)
    swap_info = rotation_analyzer.detect_swap_from_original(original_material)

    assert swap_info.is_swapped == is_swapped
    assert_two_entities_deep_almost_equal(swap_info.new_lattice, another_material.lattice)
