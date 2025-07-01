import pytest
from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.interface.grain_boundary import GrainBoundaryAnalyzer
from mat3ra.made.tools.analyze.interface.utils.holders import MatchedSubstrateFilmConfigurationHolder
from mat3ra.made.tools.build.grain_boundary.builders import GrainBoundaryBuilder
from mat3ra.made.tools.build.grain_boundary.configuration import (
    GrainBoundaryConfiguration,
)
from mat3ra.made.tools.build.grain_boundary.helpers import (
    create_grain_boundary_planar,
)
from mat3ra.made.tools.build.slab.configurations import SlabConfiguration
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from .fixtures.bulk import BULK_Si_CONVENTIONAL, BULK_Si_PRIMITIVE


@pytest.mark.parametrize(
    "phase_1_material, phase_2_material, phase_1_miller, phase_2_miller, expected_elements_min",
    [
        (
            BULK_Si_PRIMITIVE,
            BULK_Si_PRIMITIVE,
            (0, 0, 1),
            (0, 0, 1),
            8,  # Minimum expected elements
        ),
    ],
)
def test_grain_boundary_analyzer(
    phase_1_material, phase_2_material, phase_1_miller, phase_2_miller, expected_elements_min
):
    """Test the GrainBoundaryAnalyzer functionality."""
    phase_1 = Material.create(phase_1_material)
    phase_2 = Material.create(phase_2_material)

    analyzer = GrainBoundaryAnalyzer(
        phase_1_material=phase_1,
        phase_2_material=phase_2,
        phase_1_miller_indices=phase_1_miller,
        phase_2_miller_indices=phase_2_miller,
        phase_1_thickness=1,
        phase_2_thickness=1,
    )

    # Test that we get match holders
    match_holders = analyzer.grain_boundary_match_holders
    assert len(match_holders) > 0

    # Test that we can get configurations
    configs = analyzer.get_grain_boundary_configurations()
    assert len(configs) > 0

    # Test getting configuration by match ID
    config = analyzer.get_grain_boundary_configuration_by_match_id(0)
    assert isinstance(config, MatchedSubstrateFilmConfigurationHolder)
    assert config.match_id == 0


@pytest.mark.parametrize(
    "material_config, phase_1_miller, phase_2_miller, gap, translation_vector",
    [
        (
            BULK_Si_CONVENTIONAL,
            (0, 0, 1),
            (0, 1, 1),
            2.0,
            [0.0, 0.0],
        ),
    ],
)
def test_create_grain_boundary_planar(material_config, phase_1_miller, phase_2_miller, gap, translation_vector):
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
        match_id=0,
        max_area=150,
    )

    assert isinstance(grain_boundary, Material)
    assert len(grain_boundary.basis.elements.values) > 0
    assert "Grain Boundary" in grain_boundary.name


def test_grain_boundary_builder():
    """Test the GrainBoundaryBuilder directly."""
    max_area = 150
    phase_1_material = Material.create(BULK_Si_CONVENTIONAL)
    phase_2_material = Material.create(BULK_Si_CONVENTIONAL)

    # Create analyzer and get configuration
    analyzer = GrainBoundaryAnalyzer(
        phase_1_material=phase_1_material,
        phase_2_material=phase_2_material,
        phase_1_miller_indices=(0, 0, 1),
        phase_2_miller_indices=(1, 1, 1),
        phase_1_thickness=1,
        phase_2_thickness=1,
        max_area=max_area,
    )

    strained_config = analyzer.get_grain_boundary_configuration_by_match_id(0)

    # Use from_parameters to ensure gap logic is applied
    config = GrainBoundaryConfiguration.from_parameters(
        phase_1_configuration=strained_config.substrate_configuration,
        phase_2_configuration=strained_config.film_configuration,
        xy_shift=[0.0, 0.0],
        gap=3.0,
    )

    # Build grain boundary
    builder = GrainBoundaryBuilder()
    grain_boundary = builder.get_material(config)

    assert isinstance(grain_boundary, Material)
    assert len(grain_boundary.basis.elements.values) > 0
