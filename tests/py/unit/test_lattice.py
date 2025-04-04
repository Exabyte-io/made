from mat3ra.made.lattice import Lattice, LatticeVector
from mat3ra.utils import assertion as assertion_utils

DEFAULT_UNITS = Lattice.__units_default__
DEFAULT_TYPE = Lattice.__type_default__


def test_lattice_creation():
    lattice = Lattice(a=2.0, b=3.0, c=4.0)
    assert lattice.a == 2.0
    assert lattice.b == 3.0
    assert lattice.c == 4.0
    assert lattice.alpha == 90.0
    assert lattice.beta == 90.0
    assert lattice.gamma == 90.0
    assert lattice.units == DEFAULT_UNITS
    assert lattice.type == DEFAULT_TYPE


def test_lattice_get_vectors():
    lattice = Lattice(a=2.0, b=3.0, c=4.0)
    expected_vectors = [[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]]
    assert expected_vectors == lattice.vector_arrays


def test_lattice_vectors_access():
    lattice = Lattice(a=2.0, b=3.0, c=4.0)

    # Test individual vector access
    assert isinstance(lattice.vectors.a, LatticeVector)
    assert isinstance(lattice.vectors.b, LatticeVector)
    assert isinstance(lattice.vectors.c, LatticeVector)

    # Test vector arrays access
    arrays = lattice.vector_arrays
    assert len(arrays) == 3
    assert arrays[0] == lattice.vectors.a
    assert arrays[1] == lattice.vectors.b
    assert arrays[2] == lattice.vectors.c

    # Test get value of vectors, without rounding
    assertion_utils.assert_deep_almost_equal(lattice.vectors.a.root, [2.0, 0.0, 0.0])


def test_lattice_from_vectors():
    lattice = Lattice.from_vectors_array(vectors=[[2.0, 0, 0], [0, 3.0, 0], [0, 0, 4.0]])
    assert lattice.a == 2.0
    assert lattice.b == 3.0
    assert lattice.c == 4.0
    assert lattice.alpha == 90.0
    assert lattice.beta == 90.0
    assert lattice.gamma == 90.0
    assert lattice.units == DEFAULT_UNITS
    assert lattice.type.value == DEFAULT_TYPE
    assert lattice.cell_volume == 24.0
    assert lattice.vector_arrays == [[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]]


def test_lattice_get_scaled_by_matrix():
    original_lattice = Lattice(a=2.0, b=3.0, c=4.0)
    matrix = [[1.5, 0, 0], [0, 1.5, 0], [0, 0, 0.5]]
    expected_vector_values = [[3.0, 0.0, 0.0], [0.0, 4.5, 0.0], [0.0, 0.0, 2.0]]
    lattice = original_lattice.get_scaled_by_matrix(matrix)
    assert lattice.a == 3.0
    assert lattice.b == 4.5
    assert lattice.c == 2.0
    assert lattice.alpha == 90.0
    assert lattice.beta == 90.0
    assert lattice.gamma == 90.0
    assert lattice.units == DEFAULT_UNITS
    assert lattice.type.value == DEFAULT_TYPE
    assert lattice.cell_volume == 27.0
    assertion_utils.assert_deep_almost_equal(lattice.vector_arrays, expected_vector_values)


# to test: create, calculate_vectors, from_vectors, get_lattice_type, clone
# from_vectors, to_dict, cell, cell_volume, scale_by_matrix, update_from_lattice,
