from mat3ra.made.material import Material
from mat3ra.made.tools.build.grain_boundary import (
    SlabGrainBoundaryBuilder,
    SlabGrainBoundaryConfiguration,
    create_grain_boundary,
)
from mat3ra.made.tools.build.interface import ZSLStrainMatchingInterfaceBuilderParameters
from mat3ra.made.tools.build.slab import SlabConfiguration, get_terminations
from mat3ra.utils import assertion as assertion_utils


def test_slab_grain_boundary_builder():
    material = Material(Material.default_config)
    phase_1_configuration = SlabConfiguration(
        bulk=material,
        vacuum=0,
        thickness=2,
        miller_indices=(0, 0, 1),
    )

    phase_2_configuration = SlabConfiguration(
        bulk=material,
        vacuum=0,
        thickness=2,
        miller_indices=(0, 0, 1),
    )

    termination1 = get_terminations(phase_1_configuration)[0]
    termination2 = get_terminations(phase_2_configuration)[0]

    slab_config = SlabConfiguration(
        vacuum=1,
        miller_indices=(0, 0, 1),
        thickness=2,
        xy_supercell_matrix=[[1, 0], [0, 1]],
    )

    config = SlabGrainBoundaryConfiguration(
        phase_1_configuration=phase_1_configuration,
        phase_2_configuration=phase_2_configuration,
        phase_1_termination=termination1,
        phase_2_termination=termination2,
        gap=3.0,
        slab_configuration=slab_config,
    )

    builder_params = ZSLStrainMatchingInterfaceBuilderParameters(max_area=50)
    builder = SlabGrainBoundaryBuilder(build_parameters=builder_params)
    gb = create_grain_boundary(config, builder)
    expected_lattice_vectors = [
        [25.140673461, 0.0, 0.0],
        [0.0, 3.867, 0.0],
        [0.0, 0.0, 8.734],
    ]
    expected_coordinate_15 = [0.777190818, 0.5, 0.332064346]

    assert len(gb.basis.elements.values) == 32
    assertion_utils.assert_deep_almost_equal(expected_coordinate_15, gb.basis.coordinates.values[15])
    assertion_utils.assert_deep_almost_equal(expected_lattice_vectors, gb.lattice.vector_arrays)
