import json

from mat3ra.made.basis import Basis
from mat3ra.made.lattice import Lattice
from mat3ra.made.material import Material
from mat3ra.utils import assertion as assertion_utils

from mat3ra.made.utils import ArrayWithIds

REFERENCE_OBJECT_1 = {"key1": "value1", "key2": "value2"}


def test_create():
    material = Material.create_default()
    assert isinstance(material.basis, Basis)
    assert isinstance(material.basis.coordinates, ArrayWithIds)
    assert isinstance(material.lattice, Lattice)
    assert material.name == Material.__default_config__["name"]


def test_material_to_json():
    material = Material.create_default()
    # Remove all keys that are null in the config
    material_dict = json.loads(material.to_json())
    material_dict_copy = material_dict.copy()
    for key in material_dict:
        if not key in Material.__default_config__:
            material_dict_copy.pop(key)
    actual_data = json.dumps(material_dict_copy)
    expected_data = json.dumps(Material.__default_config__)
    assertion_utils.assert_deep_almost_equal(expected_data, actual_data)


def test_basis_to_json():
    material = Material.create_default()
    basis = material.basis
    expected_basis_config = {**Material.__default_config__["basis"], "labels": []}
    assertion_utils.assert_deep_almost_equal(expected_basis_config, basis.to_json())


def test_basis_cell_lattice_sync():
    """Test synchronization between basis.cell and material.lattice"""
    material = Material.create_default()

    # Change lattice vectors
    new_vectors = [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]]
    material.set_new_lattice_vectors(*new_vectors)

    # Verify basis.cell matches new lattice vectors
    assertion_utils.assert_deep_almost_equal(new_vectors, material.basis.cell.vectors_as_array)

    assertion_utils.assert_deep_almost_equal(new_vectors, material.lattice.vector_arrays)
