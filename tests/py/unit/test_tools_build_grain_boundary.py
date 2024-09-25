from mat3ra.made.material import Material
from mat3ra.made.tools.build.grain_boundary import (
    SurfaceGrainBoundaryBuilder,
    SurfaceGrainBoundaryBuilderParameters,
    SurfaceGrainBoundaryConfiguration,
)
from mat3ra.utils import assertion as assertion_utils

from .fixtures import GRAPHENE


def test_create_surface_grain_boundary():
    config = SurfaceGrainBoundaryConfiguration(
        film=Material(GRAPHENE),
        twist_angle=13.0,
        gap=2.0,
    )

    builder_params = SurfaceGrainBoundaryBuilderParameters(
        max_repetition_int=5,
        angle_tolerance=0.5,
        return_first_match=True,
        distance_tolerance=1.0,
    )
    builder = SurfaceGrainBoundaryBuilder(build_parameters=builder_params)

    gb = builder.get_materials(config)

    expected_cell_vectors = [
        [23.509344266, 0.0, 0.0],
        [5.377336066500001, 9.313819276550575, 0.0],
        [0.0, 0.0, 20.0],
    ]

    assert len(gb) == 1
    assertion_utils.assert_deep_almost_equal(expected_cell_vectors, gb[0].basis.cell.vectors_as_array)
