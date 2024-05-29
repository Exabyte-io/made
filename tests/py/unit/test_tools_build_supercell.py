from mat3ra.made.material import Material
from mat3ra.made.tools.build.supercell import create_supercell
from mat3ra.utils import assertion as assertion_utils

si_supercell_2x2 = Material.create(
    {
        "lattice": {"a": 5.131, "b": 5.131, "c": 5.131, "alpha": 90, "beta": 90, "gamma": 90},
        "basis": {
            "coordinates": [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]],
            "elements": ["Si", "Si"],
            "units": "crystal",
        },
    }
)


def test_create_supercell():
    material = Material.create(Material.default_config)
    supercell_matrix = [[2, 0], [0, 2]]
    supercell_material = create_supercell(material, supercell_matrix)

    assertion_utils.assert_deep_almost_equal(supercell_material, si_supercell_2x2)
