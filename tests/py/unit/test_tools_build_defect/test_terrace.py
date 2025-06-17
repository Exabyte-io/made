import sys

import pytest
from mat3ra.made.lattice import COORDINATE_TOLERANCE
from mat3ra.made.material import Material
from mat3ra.made.tools.build.defect.builders import TerraceSlabDefectBuilder
from mat3ra.made.tools.build.defect.configuration import TerraceSlabDefectConfiguration
from mat3ra.utils import assertion as assertion_utils
from unit.fixtures.slab import SI_CONVENTIONAL_SLAB_001


@pytest.mark.skip(reason="we'll fix before epic-7623 is merged")
@pytest.mark.parametrize(
    "crystal_config, cut_direction, pivot_coordinate, num_added_layers, expected_coords_platform",
    [
        (
            SI_CONVENTIONAL_SLAB_001,
            [1, 0, 0],
            [0.5, 0.5, 0.5],
            1,
            {"darwin": [0.627786405, 0.75, 0.671264194], "other": [0.591583068, 0.75, 0.716534426]},
        )
    ],
)
def test_create_terrace(crystal_config, cut_direction, pivot_coordinate, num_added_layers, expected_coords_platform):
    crystal = Material.create(crystal_config)
    config = TerraceSlabDefectConfiguration(
        crystal=crystal,
        cut_direction=cut_direction,
        pivot_coordinate=pivot_coordinate,
        number_of_added_layers=num_added_layers,
    )
    new_slab = TerraceSlabDefectBuilder().get_material(configuration=config)
    if sys.platform == "darwin":
        coordinate_expected = expected_coords_platform["darwin"]
    else:
        coordinate_expected = expected_coords_platform["other"]
    defect_coordinate = new_slab.basis.coordinates.values[-1]  # Use last atom (index 59) instead of previous index 35
    atol = 10 ** (-COORDINATE_TOLERANCE)
    assertion_utils.assert_deep_almost_equal(coordinate_expected, defect_coordinate, atol=atol) 