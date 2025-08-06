import numpy as np
import pytest
from mat3ra.made.material import Material
from mat3ra.made.tools.build_components.entities.reusable.three_dimensional.supercell.helpers import create_supercell
from unit.fixtures.supercell import SI_SUPERCELL_2X2X1

from .fixtures.bulk import BULK_Si_PRIMITIVE
from .utils import assert_two_entities_deep_almost_equal


@pytest.mark.parametrize(
    "material_config, supercell_matrix, expected_material_config",
    [
        (
            BULK_Si_PRIMITIVE,
            [[2, 0, 0], [0, 2, 0], [0, 0, 1]],
            SI_SUPERCELL_2X2X1,
        ),
    ],
)
def test_create_supercell(material_config, supercell_matrix, expected_material_config):
    material = Material.create(material_config)
    supercell_material = create_supercell(material, supercell_matrix)

    assert_two_entities_deep_almost_equal(supercell_material, expected_material_config)

    expected_vectors = np.array(supercell_matrix) @ np.array(material.lattice.vector_arrays)
    actual_vectors = np.array(supercell_material.lattice.vector_arrays)
    assert np.allclose(actual_vectors, expected_vectors, atol=1e-6), "Supercell vectors do not match expected values."
