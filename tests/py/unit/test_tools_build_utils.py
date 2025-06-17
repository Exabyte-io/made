import math

import pytest
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.made.material import Material
from mat3ra.made.tools.operations.core.binary import stack_two_materials

from .fixtures.bulk import BULK_Si_PRIMITIVE


@pytest.mark.parametrize(
    "material1_config, material2_config, stacking_axis, expected_a, expected_b, expected_c",
    [
        (BULK_Si_PRIMITIVE, BULK_Si_PRIMITIVE, AxisEnum.z, 3.867, 3.867, 7.734),
        (BULK_Si_PRIMITIVE, BULK_Si_PRIMITIVE, AxisEnum.x, 7.734, 3.867, 3.867),
    ],
)
def test_stack_two_materials(
    material1_config: dict,
    material2_config: dict,
    stacking_axis: AxisEnum,
    expected_a: float,
    expected_b: float,
    expected_c: float,
):
    material1 = Material.create(material1_config)
    material2 = Material.create(material2_config)

    original_atom_count1 = len(material1.basis.coordinates.values)
    original_atom_count2 = len(material2.basis.coordinates.values)

    stacked_material = stack_two_materials(material1, material2, stacking_axis)

    assert len(stacked_material.basis.coordinates.values) == original_atom_count1 + original_atom_count2

    assert math.isclose(stacked_material.lattice.a, expected_a)
    assert math.isclose(stacked_material.lattice.b, expected_b)
    assert math.isclose(stacked_material.lattice.c, expected_c)
