import math
import pytest
from unittest.mock import MagicMock

from mat3ra.made.cell.primitive_cell import get_primitive_lattice_vectors_from_config
from tests.py.unit.utils import assert_two_entities_deep_almost_equal


# Helper to mock LatticeSchema since it relies on mat3ra.esse
def create_mock_lattice(l_type=None, a=1.0, b=1.0, c=1.0, alpha=90.0, beta=90.0, gamma=90.0):
    lattice = MagicMock()
    if l_type is not None:
        lattice.type.value = l_type
    else:
        lattice.type = None

    lattice.a = a
    lattice.b = b
    lattice.c = c
    lattice.alpha = alpha
    lattice.beta = beta
    lattice.gamma = gamma
    return lattice


class TestPrimitiveCell:

    def test_primitive_lattice_vectors_na4cl4(self):
        """
        Replicates the TypeScript test for Na4Cl4 (CUB cell).
        Translates from: tests/js/cell/primitive_cell.ts
        """
        lattice = create_mock_lattice(l_type="CUB", a=5.691694)

        actual_primitive_cell = get_primitive_lattice_vectors_from_config(lattice)
        expected_primitive_cell = [
            [5.691694, 0.0, 0.0],
            [0.0, 5.691694, 0.0],
            [0.0, 0.0, 5.691694],
        ]

        assert_two_entities_deep_almost_equal(actual_primitive_cell, expected_primitive_cell)

    @pytest.mark.parametrize(
        "l_type, kwargs, expected",
        [
            # Cubic
            ("CUB", {"a": 2.0}, [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]]),
            # Face-Centered Cubic
            ("FCC", {"a": 2.0}, [[0.0, 1.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 0.0]]),
            # Body-Centered Cubic
            ("BCC", {"a": 2.0}, [[-1.0, 1.0, 1.0], [1.0, -1.0, 1.0], [1.0, 1.0, -1.0]]),
            # Tetragonal
            ("TET", {"a": 2.0, "c": 3.0}, [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]]),
            # Body-Centered Tetragonal
            ("BCT", {"a": 2.0, "c": 3.0}, [[-1.0, 1.0, 1.5], [1.0, -1.0, 1.5], [1.0, 1.0, -1.5]]),
            # Orthorhombic
            ("ORC", {"a": 2.0, "b": 3.0, "c": 4.0}, [[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]]),
            # Face-Centered Orthorhombic
            ("ORCF", {"a": 2.0, "b": 3.0, "c": 4.0}, [[0.0, 1.5, 2.0], [1.0, 0.0, 2.0], [1.0, 1.5, 0.0]]),
            # Body-Centered Orthorhombic
            ("ORCI", {"a": 2.0, "b": 3.0, "c": 4.0}, [[-1.0, 1.5, 2.0], [1.0, -1.5, 2.0], [1.0, 1.5, -2.0]]),
            # Base-Centered Orthorhombic
            ("ORCC", {"a": 2.0, "b": 3.0, "c": 4.0}, [[1.0, 1.5, 0.0], [-1.0, 1.5, 0.0], [0.0, 0.0, 4.0]]),
            # Hexagonal
            ("HEX", {"a": 2.0, "c": 3.0}, [[1.0, -math.sqrt(3), 0.0], [1.0, math.sqrt(3), 0.0], [0.0, 0.0, 3.0]]),
            # Rhombohedral (example angles)
            (
                "RHL",
                {"a": 2.0, "alpha": 60.0},
                [[1.7320508, -1.0, 0.0], [1.7320508, 1.0, 0.0], [1.1547005, 0.0, 1.6329931]],
            ),
            # Monoclinic
            (
                "MCL",
                {"a": 2.0, "b": 3.0, "c": 4.0, "alpha": 60.0},
                [[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 2.0, 3.4641016]],
            ),
            # Base-Centered Monoclinic
            (
                "MCLC",
                {"a": 2.0, "b": 3.0, "c": 4.0, "alpha": 60.0},
                [[1.0, 1.5, 0.0], [-1.0, 1.5, 0.0], [0.0, 2.0, 3.4641016]],
            ),
            # Triclinic
            (
                "TRI",
                {"a": 2.0, "b": 3.0, "c": 4.0, "alpha": 90.0, "beta": 90.0, "gamma": 90.0},
                [[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]],
            ),
        ],
    )
    def test_primitive_lattice_generators(self, l_type, kwargs, expected):
        """
        Validates the math generators for all 14 Bravais lattice scenarios.
        """
        lattice = create_mock_lattice(l_type=l_type, **kwargs)
        result = get_primitive_lattice_vectors_from_config(lattice)
        assert_two_entities_deep_almost_equal(result, expected)

    def test_fallback_to_tri_when_type_is_none(self):
        """
        Tests the behavior when lattice_config.type is None.
        Should default to the 'TRI' logic.
        """
        # A pseudo-cubic configuration passed as TRI
        lattice = create_mock_lattice(l_type=None, a=2.0, b=2.0, c=2.0, alpha=90.0, beta=90.0, gamma=90.0)
        result = get_primitive_lattice_vectors_from_config(lattice)
        expected = [
            [2.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 2.0],
        ]
        assert_two_entities_deep_almost_equal(result, expected)

    def test_unsupported_lattice_type_raises_value_error(self):
        """
        Ensures a ValueError is raised for an unknown lattice type mapping.
        """
        lattice = create_mock_lattice(l_type="INVALID_TYPE")

        with pytest.raises(ValueError) as exc_info:
            get_primitive_lattice_vectors_from_config(lattice)

        assert "Unsupported lattice type 'INVALID_TYPE'" in str(exc_info.value)
