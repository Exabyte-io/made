import pytest
from mat3ra.made.material import Material
from mat3ra.made.tools.build.supercell import create_supercell
from unit.fixtures.supercell import SI_SUPERCELL_2X2X1

from .fixtures.bulk import BULK_Si_PRIMITIVE, BULK_Si_CONVENTIONAL
from .fixtures.slab import SLAB_SrTiO3_011_TERMINATION_O2
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
