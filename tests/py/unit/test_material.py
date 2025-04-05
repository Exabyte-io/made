import numpy as np
from mat3ra.made.basis import Basis, Coordinates
from mat3ra.made.lattice import Lattice
from mat3ra.made.material import Material
from mat3ra.utils import assertion as assertion_utils
from unit.fixtures.cell import SI_CONVENTIONAL_CELL
from unit.utils import assert_two_entities_deep_almost_equal


def test_create_default():
    material = Material.create_default()
    assert isinstance(material.basis, Basis)
    assert isinstance(material.basis.coordinates, Coordinates)
    assert isinstance(material.lattice, Lattice)
    assert material.name == Material.__default_config__["name"]


def test_create():
    material = Material.create(SI_CONVENTIONAL_CELL)
    assert isinstance(material.basis, Basis)
    assert isinstance(material.basis.coordinates, Coordinates)
    assert isinstance(material.lattice, Lattice)
    assert_two_entities_deep_almost_equal(material, SI_CONVENTIONAL_CELL)


def test_create_with_cell_as_list():
    # The key cell should be ignored and Basis.Cell created from Lattice by Material
    cell = [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]
    config = {**Material.__default_config__, "basis": {**Material.__default_config__["basis"], "cell": cell}}

    material = Material.create(config)
    assert isinstance(material.basis, Basis)
    assert material.basis.cell.vector_arrays == material.lattice.vector_arrays


def test_material_to_json():
    material = Material.create_default()
    # Remove all keys that are null in the config
    assert_two_entities_deep_almost_equal(material, Material.__default_config__)


def test_material_clone():
    material = Material.create_default()
    clone = material.clone()
    assert_two_entities_deep_almost_equal(material, clone)


def test_material_to_from_cartesian():
    material = Material.create_default()
    assert material.basis.is_in_crystal_units is True
    assert material.basis.is_in_cartesian_units is False
    material.to_cartesian()
    assert material.basis.is_in_crystal_units is False
    assert material.basis.is_in_cartesian_units is True
    assert np.allclose(material.basis.coordinates.values[1], [1.1163, 0.7893, 1.9335], atol=1e-4)
    material.to_crystal()
    assert material.basis.is_in_crystal_units is True
    assert material.basis.is_in_cartesian_units is False
    assert np.allclose(material.basis.coordinates.values[1], [0.25, 0.25, 0.25], atol=1e-4)
    assert_two_entities_deep_almost_equal(material, Material.__default_config__)


def test_material_to_cartesian_to_crystal_repeated():
    material = Material.create_default()
    material.to_cartesian()
    material.to_cartesian()
    assert np.allclose(material.basis.coordinates.values[1], [1.1163, 0.7893, 1.9335], atol=1e-4)
    material.to_crystal()
    material.to_crystal()
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
    new_lattice = Lattice.from_vectors_array(new_vectors)
    material.set_lattice(new_lattice)
    # Verify basis.cell matches new lattice vectors
    assertion_utils.assert_deep_almost_equal(new_vectors, material.basis.cell.vector_arrays)
    assertion_utils.assert_deep_almost_equal(new_vectors, material.lattice.vector_arrays)
    # Verify basis coordinates are still correct
