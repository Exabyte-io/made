from mat3ra.made.basis import Basis
from mat3ra.made.lattice import Lattice
from mat3ra.made.material import Material
from mat3ra.utils import assertion as assertion_utils

REFERENCE_OBJECT_1 = {"key1": "value1", "key2": "value2"}


def test_create():
    material = Material.create(Material.default_config)
    assert isinstance(material.basis, Basis)
    assert isinstance(material.lattice, Lattice)
    assert material.name == Material.default_config["name"]


def test_material_to_json():
    material = Material.create(Material.default_config)
    assertion_utils.assert_deep_almost_equal(Material.default_config, material.to_json())


def test_basis_to_json():
    material = Material.create(Material.default_config)
    basis = material.basis
    expected_basis_config = {**Material.default_config["basis"], "labels": []}
    assertion_utils.assert_deep_almost_equal(expected_basis_config, basis.to_json())


def test_basis_cell_lattice_sync():
    """Test synchronization between basis.cell and material.lattice"""
    material = Material.create(Material.default_config)

    # Change lattice vectors
    new_vectors = [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]]
    material.set_new_lattice_vectors(*new_vectors)

    # Verify basis.cell matches new lattice vectors
    assertion_utils.assert_deep_almost_equal(new_vectors, material.basis.cell.vectors_as_array)

    assertion_utils.assert_deep_almost_equal(new_vectors, material.lattice.vector_arrays)
