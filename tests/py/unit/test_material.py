import json

from mat3ra.made.basis import Basis
from mat3ra.made.lattice import Lattice
from mat3ra.made.material import Material
from mat3ra.made.utils import ArrayWithIds
from mat3ra.utils import assertion as assertion_utils

REFERENCE_OBJECT_1 = {"key1": "value1", "key2": "value2"}


def assert_two_entities_deep_almost_equal(entity1, entity2):
    dict_1 = json.loads(entity1.to_json())
    dict_1_copy = dict_1.copy()
    for key in dict_1:
        if key not in entity2:
            dict_1_copy.pop(key)
    actual_data = json.dumps(dict_1_copy)
    expected_data = json.dumps(entity2)
    assertion_utils.assert_deep_almost_equal(expected_data, actual_data)


def test_create():
    material = Material.create_default()
    assert isinstance(material.basis, Basis)
    assert isinstance(material.basis.coordinates, ArrayWithIds)
    assert isinstance(material.lattice, Lattice)
    assert material.name == Material.__default_config__["name"]


def test_material_to_json():
    material = Material.create_default()
    # Remove all keys that are null in the config
    assert_two_entities_deep_almost_equal(material, Material.__default_config__)


def test_basis_to_json():
    material = Material.create_default()
    basis = material.basis
    assert_two_entities_deep_almost_equal(basis, Material.__default_config__["basis"])


def test_basis_cell_lattice_sync():
    """Test synchronization between basis.cell and material.lattice"""
    material = Material.create_default()
    # Change lattice vectors
    new_vectors = [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]]
    material.set_new_lattice_vectors(*new_vectors)
    # Verify basis.cell matches new lattice vectors
    assertion_utils.assert_deep_almost_equal(new_vectors, material.basis.cell.vectors_as_array)
    assertion_utils.assert_deep_almost_equal(new_vectors, material.lattice.vector_arrays)
