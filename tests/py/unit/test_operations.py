import math

import numpy as np
import pytest
from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.made.material import Material
from mat3ra.made.tools.build_components.operations.core.modifications.perturb import SineWavePerturbationFunctionHolder
from mat3ra.made.tools.operations.core.binary import stack_two_materials
from mat3ra.made.tools.operations.core.unary import perturb, strain
from unit.fixtures.bulk import BULK_Si_CONVENTIONAL
from unit.fixtures.strain import BULK_Si_CONVENTIONAL_STRAINED
from unit.utils import assert_two_entities_deep_almost_equal

from .fixtures.bulk import BULK_Si_PRIMITIVE
from .fixtures.nanoribbon.nanoribbon import GRAPHENE_ZIGZAG_NANORIBBON

STRAIN_TEST_CASES = [
    (BULK_Si_CONVENTIONAL, [[1.1, 0, 0], [0, 1.1, 0], [0, 0, 1.0]], BULK_Si_CONVENTIONAL_STRAINED),
]


@pytest.mark.parametrize("material_config, strain_matrix, expected_material_config", STRAIN_TEST_CASES)
def test_strain(material_config, strain_matrix, expected_material_config):
    material = Material.create(material_config)
    strain_matrix_obj = Matrix3x3Schema(root=strain_matrix)
    expected_material = Material.create(expected_material_config)

    strained_material = strain(material, strain_matrix_obj)

    assert_two_entities_deep_almost_equal(strained_material, expected_material)


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


@pytest.mark.parametrize(
    "material_config, perturbation_function, expected_coord_changes",
    [
        (
            GRAPHENE_ZIGZAG_NANORIBBON,
            SineWavePerturbationFunctionHolder(amplitude=0.1, wavelength=1.0, phase=0.0, axis="x"),
            [0.0, 0.0, -0.0383],
        ),
    ],
)
def test_perturb(material_config, perturbation_function, expected_coord_changes):
    material = Material.create(material_config)
    original_coords = [coord[:] for coord in material.basis.coordinates.values]

    perturbed_material = perturb(material, perturbation_function)
    actual_delta = np.array(original_coords[0]) - np.array(perturbed_material.basis.coordinates.values[0])
    assert np.isclose(actual_delta, expected_coord_changes, atol=1e-4).all()
