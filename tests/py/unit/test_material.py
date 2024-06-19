import json

from mat3ra.made.basis.basis import Basis
from mat3ra.made.lattice.lattice import Lattice
from mat3ra.made.material import Material
from mat3ra.utils import assertion as assertion_utils

REFERENCE_OBJECT_1 = {"key1": "value1", "key2": "value2"}
material = Material.create(Material.default_config)


def test_create():
    assert isinstance(material.basis, Basis)
    assert isinstance(material.lattice, Lattice)
    assert material.to_json() == Material.default_config
    assert material.name == Material.default_config["name"]


def test_material_to_json():
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
    basis = material.basis
    expected_basis_config = {**Material.default_config["basis"], "cell": None, "labels": []}
    expected_basis_json = json.loads(json.dumps(expected_basis_config))
    assertion_utils.assert_deep_almost_equal(expected_basis_json, basis.to_json())
