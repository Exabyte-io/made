from mat3ra.made.material import Material
from mat3ra.made.tools.build.supercell import create_supercell
from unit.fixtures.supercell import SI_SUPERCELL_2X2X1

from .utils import assert_two_entities_deep_almost_equal


def test_create_supercell():
    material = Material.create_default()
    supercell_matrix = [[2, 0, 0], [0, 2, 0], [0, 0, 1]]
    supercell_material = create_supercell(material, supercell_matrix)

    assert_two_entities_deep_almost_equal(supercell_material, SI_SUPERCELL_2X2X1)
