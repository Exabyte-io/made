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


def test_create_empty():
    """Test that creating an empty material results in empty basis"""
    material = Material.create_empty()
    assert material.basis.elements.values == []
    assert material.basis.coordinates.values == []


def test_create_empty_default_params():
    """Test default parameters when creating empty material"""
    material = Material.create_empty()
    assert material.name == "New Material"
    assert material.lattice is not None
    assert material.basis is not None


def test_create_empty_custom_lattice():
    """Test custom lattice parameters when creating empty material"""
    material = Material.create_empty(
        a=2.0, b=3.0, c=4.0, alpha=80.0, beta=85.0, gamma=95.0, lattice_type="TRI"
    )
    assert material.lattice.type == "TRI"
    assert material.lattice.a == 2.0
    assert material.lattice.b == 3.0
    assert material.lattice.c == 4.0
    assert material.lattice.alpha == 80.0
    assert material.lattice.beta == 85.0
    assert material.lattice.gamma == 95.0


def test_lattice_vectors_access():
    lattice = Lattice(a=2.0, b=3.0, c=4.0)

    # Test individual vector access
    assert isinstance(lattice.vectors.a, list)
    assert isinstance(lattice.vectors.b, list)
    assert isinstance(lattice.vectors.c, list)

    # Test vector arrays access
    arrays = lattice.vector_arrays
    assert len(arrays) == 3
    assert arrays[0] == lattice.vectors.a
    assert arrays[1] == lattice.vectors.b
    assert arrays[2] == lattice.vectors.c
