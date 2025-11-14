import pytest

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.lattice import LatticeMaterialAnalyzer
from mat3ra.made.tools.analyze.lattice_swap_analyzer import MaterialLatticeSwapAnalyzer
from .fixtures.interface.gaas_dia import (
    GALLIUM_ARSENIDE_DIAMOND_INTERFACE,
    GALLIUM_ARSENIDE_DIAMOND_INTERFACE_PRIMITIVE,
    GALLIUM_ARSENIDE_DIAMOND_INTERFACE_PRIMITIVE_GH_WF,
)
from .fixtures.interface.gr_ni_111_top_hcp import GRAPHENE_NICKEL_INTERFACE_TOP_HCP_GH_WF
from .fixtures.slab import (
    SLAB_SrTiO3_011_TERMINATION_O2,
    SLAB_SrTiO3_011_TERMINATION_SrTiO,
)
from .utils import assert_two_entities_deep_almost_equal, get_platform_specific_value, OSPlatform


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


def test_swap_detection():
    """
    We test that the lattice swap detection correctly identifies the transformation needed to convert
    primitive cell that is returned with swapped lattice vectors (rotated) to the correct orientation.
    Size of original and primitive materials is different, so fingerprint comparison is used.
    """
    # Case 1: Primitive material here will be "rotated" 90 degrees around X axis compared to original
    original_material = Material.create(GALLIUM_ARSENIDE_DIAMOND_INTERFACE)
    primitive_material = LatticeMaterialAnalyzer(material=original_material).material_with_primitive_lattice

    analyzer = MaterialLatticeSwapAnalyzer(material=primitive_material)
    swap_info = analyzer.detect_swap_from_original(original_material)
    assert swap_info.is_swapped is True

    corrected_primitive_material = analyzer.correct_material_to_match_target(original_material)

    expected_primitive = get_platform_specific_value(
        {
            OSPlatform.DARWIN: GALLIUM_ARSENIDE_DIAMOND_INTERFACE_PRIMITIVE,
            OSPlatform.OTHER: GALLIUM_ARSENIDE_DIAMOND_INTERFACE_PRIMITIVE_GH_WF,
        }
    )
    assert_two_entities_deep_almost_equal(corrected_primitive_material, expected_primitive)

    # Case 2: Primitive material "rotated" 180 degrees around X axis, or mirrored along z, compared to original
    original_material = Material.create(GRAPHENE_NICKEL_INTERFACE_TOP_HCP_GH_WF)
    primitive_material = LatticeMaterialAnalyzer(material=original_material).material_with_primitive_lattice

    analyzer = MaterialLatticeSwapAnalyzer(material=primitive_material)
    swap_info = analyzer.detect_swap_from_original(original_material)

    original_material.basis.set_labels_from_list([])

    corrected_primitive_material = analyzer.correct_material_to_match_target(original_material)
    print("GH MATERIAL PRIMITIVE", primitive_material)
    print("GH MATERIAL CORRECTED", corrected_primitive_material)

    assert swap_info.is_swapped is True
    assert_two_entities_deep_almost_equal(corrected_primitive_material.basis, original_material.basis)
    assert_two_entities_deep_almost_equal(corrected_primitive_material.lattice, original_material.lattice)
