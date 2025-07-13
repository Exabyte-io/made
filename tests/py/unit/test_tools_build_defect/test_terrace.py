import sys

import numpy as np
import pytest
from pymatgen.cli.pmg_analyze import analyze

from mat3ra.made.lattice import COORDINATE_TOLERANCE
from mat3ra.made.tools.build import MaterialWithBuildMetadata
from mat3ra.made.tools.build.defect.terrace.helpers import create_terrace
from mat3ra.utils import assertion as assertion_utils

from mat3ra.made.tools.modify import rotate
from mat3ra.made.tools.operations.core.unary import edit_cell, supercell
from unit.fixtures.slab import SI_CONVENTIONAL_SLAB_001


@pytest.mark.parametrize(
    "crystal_config, cut_direction, pivot_coordinate, num_added_layers, expected_coordinates_platform",
    [
        (
            SI_CONVENTIONAL_SLAB_001,
            [0, 1, 0],
            [0.5, 0.5, 0.5],
            0.25,
            {"darwin": [0.627786405, 0.75, 0.671264194], "other": [0.591583068, 0.75, 0.716534426]},
        )
    ],
)
def test_create_terrace(
    crystal_config, cut_direction, pivot_coordinate, num_added_layers, expected_coordinates_platform
):
    crystal = MaterialWithBuildMetadata.create(crystal_config)
    crystal = supercell(crystal, [[3, 0, 0], [0, 3, 0], [0, 0, 1]])
    terrace = create_terrace(
        slab=crystal,
        cut_direction=cut_direction,
        pivot_coordinate=pivot_coordinate,
        number_of_added_layers=num_added_layers,
    )
    if sys.platform == "darwin":
        coordinate_expected = expected_coordinates_platform["darwin"]
    else:
        coordinate_expected = expected_coordinates_platform["other"]
    defect_coordinate = terrace.basis.coordinates.values[-1]  # Use last atom (index 59) instead of previous index 35
    atol = 10 ** (-COORDINATE_TOLERANCE)

    assertion_utils.assert_deep_almost_equal(coordinate_expected, defect_coordinate, atol=atol)
