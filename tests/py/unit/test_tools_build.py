import pytest
from mat3ra.made.material import Material
from mat3ra.made.tools.operations.core.binary import merge_materials
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


@pytest.mark.parametrize(
    "materials_to_merge, expected_basis, expected_basis_reverse",
    [
        (
            [section, cavity],
            MERGED_CAVITY_SECTION_BASIS,
            MERGED_SECTION_CAVITY_BASIS,
        ),
        (
            [section_with_extra_atom, cavity],
            MERGED_CAVITY_SECTION_BASIS,
            MERGED_SECTION_CAVITY_BASIS,
        ),
    ],
)
def test_merge_materials(materials_to_merge, expected_basis, expected_basis_reverse):
    merged_material = merge_materials(materials_to_merge)
    merged_material_reverse = merge_materials(materials_to_merge[::-1])
    assertion_utils.assert_deep_almost_equal(merged_material.basis, expected_basis)
    assertion_utils.assert_deep_almost_equal(merged_material_reverse.basis, expected_basis_reverse)
