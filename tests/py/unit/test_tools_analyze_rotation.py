"""
Tests for MaterialRotationAnalyzer class.
"""

import numpy as np
import pytest

from mat3ra.made.lattice import Lattice
from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.lattice import LatticeMaterialAnalyzer
from mat3ra.made.tools.analyze.lattice_swap_analyzer import MaterialLatticeSwapAnalyzer
from mat3ra.made.tools.modify import wrap_to_unit_cell
from mat3ra.made.tools.operations.core.unary import rotate, translate
from .fixtures.interface.gaas_dia import GALLIUM_ARSENIDE_DIAMOND_INTERFACE
from .fixtures.interface.zsl import DIAMOND_GaAs_INTERFACE
from .fixtures.slab import (
    SLAB_SrTiO3_011_TERMINATION_O2,
    SLAB_SrTiO3_011_TERMINATION_SrTiO,
)
from .test_tools_build_interface import Si_Ge_SIMPLE_INTERFACE_TEST_CASE
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


def test_lattice_swap_detection_primitive():
    """Test that lattice swap detection correctly identifies vector permutations."""
    original_material = Material.create(GALLIUM_ARSENIDE_DIAMOND_INTERFACE)
    analyzer = LatticeMaterialAnalyzer(material=original_material)

    corrected_primitive_material = analyzer.get_material_with_primitive_lattice_standard(keep_orientation=True)
    # The corrected version should have the same lattice.c as the original
    assert abs(original_material.lattice.c - corrected_primitive_material.lattice.c) < 0.01
    assert abs(corrected_primitive_material.lattice.a - 5.63) < 0.01
    assert abs(corrected_primitive_material.lattice.b - 5.63) < 0.01


# @pytest.mark.parametrize(
#     "original_material_config, material_with_swapped_lattice_config, expected_permutation",
#     [(DIAMOND_GaAs_INTERFACE, DIAMOND_GaAs_INTERFACE, [(2, 1), (0, 1), (1, -1)])],
# )
# def test_lattice_swap_detection_between_materials(
#     original_material_config, material_with_swapped_lattice_config, expected_permutation
# ):
#     """Test lattice swap detection between materials with known permutations."""
#     # TODO: Generalize. Temporarily manually setting rotation to check for
#     original_material = Material.create(original_material_config)
#     material_with_swapped_lattice = original_material.clone()
#     original_lattice_vectors = np.array(original_material.lattice.vector_arrays)
#     swapped_lattice_vectors = [
#         original_lattice_vectors[2],
#         original_lattice_vectors[0],
#         -original_lattice_vectors[1],
#     ]
#     swapped_lattice = Lattice.from_vectors_array(swapped_lattice_vectors)
#     material_with_swapped_lattice.set_lattice(swapped_lattice)
#     material_with_swapped_lattice = rotate(material_with_swapped_lattice, axis=[0, 1, 0], angle=-90)
#
#     rotation_analyzer = MaterialLatticeSwapAnalyzer(material=material_with_swapped_lattice)
#     swap_info = rotation_analyzer.detect_swap_from_original(original_material)
#     corrected_material = rotation_analyzer.correct_material_to_match_original(original_material)
#     assert swap_info.is_swapped is True
#     assert swap_info.permutation == expected_permutation
