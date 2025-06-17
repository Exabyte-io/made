import pytest
from mat3ra.made.material import Material
from mat3ra.made.tools.build.grain_boundary import (
    SlabGrainBoundaryConfiguration,
    SurfaceGrainBoundaryBuilder,
    SurfaceGrainBoundaryBuilderParameters,
    SurfaceGrainBoundaryConfiguration,
    create_grain_boundary,
)
from mat3ra.made.tools.build.grain_boundary.builders import SlabGrainBoundaryBuilder, SlabGrainBoundaryBuilderParameters
from mat3ra.made.tools.build.slab.configuration import SlabConfiguration
from mat3ra.utils import assertion as assertion_utils

from .fixtures.bulk import BULK_Si_PRIMITIVE
from .fixtures.monolayer import GRAPHENE


@pytest.mark.skip(reason="Failing due to SlabConfiguration instantiation changes. To be fixed in epic-7623")
@pytest.mark.parametrize(
    "material_config, slab_params, gap, expected_elements_len, expected_coordinate_checks, expected_lattice_vectors",
    [
        (
            BULK_Si_PRIMITIVE,
            {"vacuum": 0, "number_of_layers": 2, "miller_indices": (0, 0, 1)},
            3.0,
            32,
            {15: [0.777190818, 0.5, 0.332064346]},
            [[25.140673461, 0.0, 0.0], [0.0, 3.867, 0.0], [0.0, 0.0, 8.734]],
        ),
    ],
)
def test_slab_grain_boundary_builder(
    material_config, slab_params, gap, expected_elements_len, expected_coordinate_checks, expected_lattice_vectors
):
    material = Material.create(material_config)
    phase_1_configuration = SlabConfiguration(bulk=material, **slab_params)
    phase_2_configuration = SlabConfiguration(bulk=material, **slab_params)

    termination1 = phase_1_configuration.get_terminations()[0]
    termination2 = phase_2_configuration.get_terminations()[0]

    slab_config_params = slab_params.copy()
    slab_config_params.update({"vacuum": 1, "xy_supercell_matrix": [[1, 0], [0, 1]]})
    slab_config = SlabConfiguration(bulk=material, **slab_config_params)

    config = SlabGrainBoundaryConfiguration(
        phase_1_configuration=phase_1_configuration,
        phase_2_configuration=phase_2_configuration,
        phase_1_termination=termination1,
        phase_2_termination=termination2,
        gap=gap,
        slab_configuration=slab_config,
    )

    builder_params = SlabGrainBoundaryBuilderParameters()
    builder = SlabGrainBoundaryBuilder(build_parameters=builder_params)
    gb = create_grain_boundary(config, builder)

    assert len(gb.basis.elements.values) == expected_elements_len
    for index, expected_coordinate in expected_coordinate_checks.items():
        assertion_utils.assert_deep_almost_equal(expected_coordinate, gb.basis.coordinates.values[index])
    assertion_utils.assert_deep_almost_equal(expected_lattice_vectors, gb.lattice.vector_arrays)


@pytest.mark.parametrize(
    "config_params, builder_params_dict, expected_cell_vectors",
    [
        (
            {"film_config": GRAPHENE, "twist_angle": 13.0, "gap": 2.0},
            {
                "max_repetition_int": 5,
                "angle_tolerance": 0.5,
                "return_first_match": True,
                "distance_tolerance": 1.0,
            },
            [
                [23.509344266, 0.0, 0.0],
                [5.377336066500001, 9.313819276550575, 0.0],
                [0.0, 0.0, 20.0],
            ],
        ),
    ],
)
def test_create_surface_grain_boundary(config_params, builder_params_dict, expected_cell_vectors):
    config_params["film"] = Material.create(config_params.pop("film_config"))
    config = SurfaceGrainBoundaryConfiguration(**config_params)
    builder_params = SurfaceGrainBoundaryBuilderParameters(**builder_params_dict)
    builder = SurfaceGrainBoundaryBuilder(build_parameters=builder_params)
    gb = builder.get_materials(config)

    assert len(gb) == 1
    assertion_utils.assert_deep_almost_equal(expected_cell_vectors, gb[0].basis.cell.vector_arrays)
