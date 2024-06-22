from mat3ra.made.basis import Basis
from mat3ra.made.lattice import Lattice
from mat3ra.made.material import Material
from mat3ra.utils import assertion as assertion_utils

REFERENCE_OBJECT_1 = {"key1": "value1", "key2": "value2"}


def test_create():
    material = Material.create(Material.default_config)
    assert isinstance(material.basis, Basis)
    assert isinstance(material.lattice, Lattice)
    assert material.to_json() == Material.default_config
    assert material.name == Material.default_config["name"]


def test_material_to_json():
    material = Material.create(Material.default_config)
    labels_array = [{"id": 0, "value": 0}, {"id": 1, "value": 1}]
    config_with_labels = {
        **Material.default_config,
        "basis": {**Material.default_config["basis"], "labels": labels_array},
    }
    expected_config = {
        **Material.default_config,
        "basis": {**Material.default_config["basis"], "cell": None, "labels": labels_array},
    }
    material.basis = Basis.from_dict(**config_with_labels["basis"])
    assertion_utils.assert_deep_almost_equal(expected_config, material.to_json())


def test_basis_to_json():
    material = Material.create(Material.default_config)
    basis = material.basis
    expected_basis_config = {**Material.default_config["basis"], "cell": None, "labels": []}
    assertion_utils.assert_deep_almost_equal(expected_basis_config, basis.to_json())
