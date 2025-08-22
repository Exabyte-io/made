import pytest

from mat3ra.made.tools.analyze.lattice import get_lattice_type
from .fixtures.bulk import BULK_Si_PRIMITIVE, BULK_Si_CONVENTIONAL
from .fixtures.interface.gr_ni_111_top_hcp import GRAPHENE_NICKEL_INTERFACE_TOP_HCP


@pytest.mark.parametrize(
    "material, expected_lattice_type",
    [
        (BULK_Si_PRIMITIVE, "FCC"),
        (BULK_Si_CONVENTIONAL, "CUB"),
        (GRAPHENE_NICKEL_INTERFACE_TOP_HCP, "HEX"),
    ],
)
def test_analyze_lattice_type(material, expected_lattice_type):
    result = get_lattice_type(material)
    assert result.value == expected_lattice_type
