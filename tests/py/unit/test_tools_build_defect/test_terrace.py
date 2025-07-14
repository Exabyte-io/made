import sys

import pytest
from mat3ra.utils import assertion as assertion_utils

from mat3ra.made.lattice import COORDINATE_TOLERANCE
from mat3ra.made.material import Material
from mat3ra.made.tools.build import MaterialWithBuildMetadata
from mat3ra.made.tools.build.defect.terrace.helpers import create_terrace
from mat3ra.made.tools.build.slab.helpers import create_slab
from mat3ra.made.tools.operations.core.unary import supercell
from unit.fixtures.bulk import BULK_Si_CONVENTIONAL
from unit.fixtures.slab import SI_CONVENTIONAL_SLAB_001
from unit.fixtures.terrace import TERRACE_SLAB_Si_001_3x3
from unit.utils import assert_two_entities_deep_almost_equal


@pytest.mark.parametrize(
    "slab_parameters, cut_direction, pivot_coordinate, num_added_layers, expected_config",
    [
        (
            {"crystal": BULK_Si_CONVENTIONAL, "number_of_layers": 2, "xy_supercell_matrix": [[3, 0], [0, 3]]},
            [0, 1, 0],
            [0.5, 0.5, 0.5],
            0.5,
            TERRACE_SLAB_Si_001_3x3,
        )
    ],
)
def test_create_terrace(
    slab_parameters,
    cut_direction,
    pivot_coordinate,
    num_added_layers,
    expected_config,
):
    slab = create_slab(
        Material.create(slab_parameters["crystal"]),
        number_of_layers=slab_parameters["number_of_layers"],
        xy_supercell_matrix=slab_parameters["xy_supercell_matrix"],
    )
    terrace = create_terrace(
        slab=slab,
        cut_direction=cut_direction,
        pivot_coordinate=pivot_coordinate,
        number_of_added_layers=num_added_layers,
    )
    # if sys.platform == "darwin":
    #     coordinate_expected = expected_coordinates_platform["darwin"]
    # else:
    #     coordinate_expected = expected_coordinates_platform["other"]
    defect_coordinate = terrace.basis.coordinates.values[-1]
    atol = 10 ** (-COORDINATE_TOLERANCE)

    assert_two_entities_deep_almost_equal(terrace, expected_config)
