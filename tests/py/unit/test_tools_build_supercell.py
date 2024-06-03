from mat3ra.made.material import Material
from mat3ra.made.tools.build.supercell import create_supercell
from mat3ra.utils import assertion as assertion_utils

from .fixtures import SI_SUPERCELL_2X2X1

si_supercell_2x2x1 = Material(SI_SUPERCELL_2X2X1)


def test_create_supercell():
    material = Material.create(Material.default_config)
    supercell_matrix = [[2, 0, 0], [0, 2, 0], [0, 0, 1]]
    supercell_material = create_supercell(material, supercell_matrix)

    assertion_utils.assert_deep_almost_equal(supercell_material.to_json(), si_supercell_2x2x1.to_json())
