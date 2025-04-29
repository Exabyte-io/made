from mat3ra.made.material import Material
from mat3ra.made.tools.build.utils import merge_materials
from mat3ra.utils import assertion as assertion_utils
from unit.fixtures.cuts import (
    CAVITY_MATERIAL_BASIS,
    FULL_MATERIAL,
    MERGED_CAVITY_SECTION_BASIS,
    MERGED_SECTION_CAVITY_BASIS,
    SECTION_MATERIAL_BASIS,
    SECTION_MATERIAL_BASIS_EXTRA_ATOM,
)

section = Material.create({**FULL_MATERIAL, **SECTION_MATERIAL_BASIS})
cavity = Material.create({**FULL_MATERIAL, **CAVITY_MATERIAL_BASIS})
section_with_extra_atom = Material.create({**FULL_MATERIAL, **SECTION_MATERIAL_BASIS_EXTRA_ATOM})


def test_merge_materials():
    merged_material = merge_materials([section, cavity])
    merged_material_reverse = merge_materials([cavity, section])
    assertion_utils.assert_deep_almost_equal(merged_material.basis, MERGED_CAVITY_SECTION_BASIS)
    assertion_utils.assert_deep_almost_equal(merged_material_reverse.basis, MERGED_SECTION_CAVITY_BASIS)


def test_resolve_close_coordinates_basis():
    merged_material = merge_materials([section, cavity])
    merged_material_reverse = merge_materials([cavity, section])
    assertion_utils.assert_deep_almost_equal(merged_material.basis, MERGED_CAVITY_SECTION_BASIS)
    assertion_utils.assert_deep_almost_equal(merged_material_reverse.basis, MERGED_SECTION_CAVITY_BASIS)


def test_resolve_close_coordinates_basis_extra_atom():
    merged_material = merge_materials([section_with_extra_atom, cavity])
    merged_material_reverse = merge_materials([cavity, section_with_extra_atom])
    assertion_utils.assert_deep_almost_equal(merged_material.basis, MERGED_CAVITY_SECTION_BASIS)
    assertion_utils.assert_deep_almost_equal(merged_material_reverse.basis, MERGED_SECTION_CAVITY_BASIS)
