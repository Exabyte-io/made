import pytest
from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.interface.grain_boundary import GrainBoundaryAnalyzer
from mat3ra.made.tools.build.grain_boundary.builders import GrainBoundaryBuilder
from mat3ra.made.tools.build.grain_boundary.configuration import GrainBoundaryConfiguration
from mat3ra.made.tools.build.grain_boundary.helpers import create_grain_boundary_planar

from .fixtures.bulk import BULK_Si_CONVENTIONAL
from .fixtures.grain_boundary import GRAIN_BOUNDARY_SI_001_011
from .utils import assert_two_entities_deep_almost_equal


@pytest.mark.parametrize(
    "material_config, phase_1_miller, phase_2_miller, gap, translation_vector, max_area, expected_material_config",
    [
        (
            BULK_Si_CONVENTIONAL,
            (0, 0, 1),
            (0, 1, 1),
            2.0,
            [2.0, 1.0],
            150,
            GRAIN_BOUNDARY_SI_001_011,
        ),
    ],
)
def test_create_grain_boundary_planar(
    material_config, phase_1_miller, phase_2_miller, gap, translation_vector, max_area, expected_material_config
):
    """Test creating a planar grain boundary."""
    phase_1 = Material.create(material_config)

    grain_boundary = create_grain_boundary_planar(
        phase_1_material=phase_1,
        phase_1_miller_indices=phase_1_miller,
        phase_2_miller_indices=phase_2_miller,
        phase_1_thickness=1,
        phase_2_thickness=1,
        translation_vector=translation_vector,
        gap=gap,
        max_area=max_area,
    )

    assert_two_entities_deep_almost_equal(grain_boundary, expected_material_config)


@pytest.mark.parametrize(
    "material_config, phase_1_miller, phase_2_miller, gap, translation_vector, max_area, expected_material_config",
    [
        (
            BULK_Si_CONVENTIONAL,
            (0, 0, 1),
            (0, 1, 1),
            2.0,
            [2.0, 1.0],
            150,
            GRAIN_BOUNDARY_SI_001_011,
        ),
    ],
)
def test_grain_boundary_builder(
    material_config, phase_1_miller, phase_2_miller, gap, translation_vector, max_area, expected_material_config
):
    phase_1_material = Material.create(material_config)
    phase_2_material = Material.create(material_config)

    # Create analyzer and get configuration
    analyzer = GrainBoundaryAnalyzer(
        phase_1_material=phase_1_material,
        phase_2_material=phase_2_material,
        phase_1_miller_indices=phase_1_miller,
        phase_2_miller_indices=phase_2_miller,
        phase_1_thickness=1,
        phase_2_thickness=1,
        max_area=max_area,
    )

    strained_config = analyzer.get_grain_boundary_configuration_by_match_id(0)

    # Use from_parameters to ensure gap logic is applied
    config = GrainBoundaryConfiguration.from_parameters(
        phase_1_configuration=strained_config.substrate_configuration,
        phase_2_configuration=strained_config.film_configuration,
        xy_shift=translation_vector,
        gap=gap,
    )

    # Build grain boundary
    builder = GrainBoundaryBuilder()
    grain_boundary = builder.get_material(config)

    assert_two_entities_deep_almost_equal(grain_boundary, expected_material_config)
