import json

import numpy as np
from mat3ra.made.basis import Basis
from mat3ra.made.lattice import Lattice
from mat3ra.made.material import Material
from mat3ra.made.utils import ArrayWithIds
from mat3ra.utils import assertion as assertion_utils

from unit.utils import assert_two_entities_deep_almost_equal

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
    assert_two_entities_deep_almost_equal(material, Material.__default_config__)


def test_material_to_from_cartesian():
    material = Material.create_default()
    material.basis.to_cartesian()
    assert np.allclose(material.basis.coordinates.values[1], [1.1163, 0.7893, 1.9335], atol=1e-4)
    material.basis.to_crystal()
    assert np.allclose(material.basis.coordinates.values[1], [0.25, 0.25, 0.25], atol=1e-4)
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
