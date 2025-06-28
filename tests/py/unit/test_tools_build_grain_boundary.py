import pytest
from mat3ra.made.material import Material
from mat3ra.made.tools.build.grain_boundary.helpers import (
    create_grain_boundary_planar,
    create_grain_boundary_planar_with_vacuum,
    create_grain_boundary,
)
from mat3ra.made.tools.build.grain_boundary.builders import (
    GrainBoundaryBuilder,
    GrainBoundaryWithVacuumBuilder,
)
from mat3ra.made.tools.build.grain_boundary.configuration import (
    GrainBoundaryConfiguration,
    GrainBoundaryWithVacuumConfiguration,
)
from mat3ra.made.tools.analyze.interface.grain_boundary import GrainBoundaryAnalyzer
from mat3ra.made.tools.analyze.interface.utils.holders import MatchedSubstrateFilmConfigurationHolder
from mat3ra.utils import assertion as assertion_utils

from .fixtures.bulk import BULK_Si_PRIMITIVE


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
    "phase_1_material, phase_2_material, phase_1_miller, phase_2_miller, gap, translation_vector",
    [
        (
            BULK_Si_PRIMITIVE,
            BULK_Si_PRIMITIVE,
            (0, 0, 1),
            (0, 0, 1),
            3.0,
            [0.0, 0.0, 0.0],
        ),
    ],
)
def test_create_grain_boundary_planar(
    phase_1_material, phase_2_material, phase_1_miller, phase_2_miller, gap, translation_vector
):
    """Test creating a planar grain boundary."""
    phase_1 = Material.create(phase_1_material)
    phase_2 = Material.create(phase_2_material)

    grain_boundary = create_grain_boundary_planar(
        phase_1_material=phase_1,
        phase_2_material=phase_2,
        phase_1_miller_indices=phase_1_miller,
        phase_2_miller_indices=phase_2_miller,
        phase_1_thickness=1,
        phase_2_thickness=1,
        translation_vector=translation_vector,
        gap=gap,
        match_id=0,
    )

    assert isinstance(grain_boundary, Material)
    assert len(grain_boundary.basis.elements.values) > 0
    assert "Grain Boundary" in grain_boundary.name


@pytest.mark.parametrize(
    "phase_1_material, phase_2_material, phase_1_miller, phase_2_miller, slab_miller, gap, vacuum",
    [
        (
            BULK_Si_PRIMITIVE,
            BULK_Si_PRIMITIVE,
            (0, 0, 1),
            (0, 0, 1),
            (0, 0, 1),
            3.0,
            10.0,
        ),
    ],
)
def test_create_grain_boundary_planar_with_vacuum(
    phase_1_material, phase_2_material, phase_1_miller, phase_2_miller, slab_miller, gap, vacuum
):
    """Test creating a grain boundary with vacuum."""
    phase_1 = Material.create(phase_1_material)
    phase_2 = Material.create(phase_2_material)

    grain_boundary_slab = create_grain_boundary_planar_with_vacuum(
        phase_1_material=phase_1,
        phase_2_material=phase_2,
        phase_1_miller_indices=phase_1_miller,
        phase_2_miller_indices=phase_2_miller,
        slab_miller_indices=slab_miller,
        phase_1_thickness=1,
        phase_2_thickness=1,
        slab_thickness=1,
        gap=gap,
        vacuum=vacuum,
        match_id=0,
    )

    assert isinstance(grain_boundary_slab, Material)
    assert len(grain_boundary_slab.basis.elements.values) > 0
    assert "Grain Boundary with Vacuum" in grain_boundary_slab.name


def test_grain_boundary_builder():
    """Test the GrainBoundaryBuilder directly."""
    phase_1 = Material.create(BULK_Si_PRIMITIVE)
    phase_2 = Material.create(BULK_Si_PRIMITIVE)

    # Create analyzer and get configuration
    analyzer = GrainBoundaryAnalyzer(
        phase_1_material=phase_1,
        phase_2_material=phase_2,
        phase_1_miller_indices=(0, 0, 1),
        phase_2_miller_indices=(0, 0, 1),
    )

    strained_config = analyzer.get_grain_boundary_configuration_by_match_id(0)

    # Create configuration
    config = GrainBoundaryConfiguration(
        phase_1_configuration=strained_config,
        phase_2_configuration=strained_config,
        translation_vector=[0.0, 0.0, 0.0],
        gap=3.0,
    )

    # Build grain boundary
    builder = GrainBoundaryBuilder()
    grain_boundary = builder.get_material(config)

    assert isinstance(grain_boundary, Material)
    assert len(grain_boundary.basis.elements.values) > 0


def test_grain_boundary_with_vacuum_builder():
    """Test the GrainBoundaryWithVacuumBuilder directly."""
    # First create a grain boundary material
    phase_1 = Material.create(BULK_Si_PRIMITIVE)
    phase_2 = Material.create(BULK_Si_PRIMITIVE)

    grain_boundary_material = create_grain_boundary_planar(
        phase_1_material=phase_1,
        phase_2_material=phase_2,
        phase_1_miller_indices=(0, 0, 1),
        phase_2_miller_indices=(0, 0, 1),
        match_id=0,
    )

    # Create configuration for grain boundary with vacuum
    config = GrainBoundaryWithVacuumConfiguration.from_grain_boundary_material(
        grain_boundary_material=grain_boundary_material,
        miller_indices=(0, 0, 1),
        number_of_layers=1,
        vacuum=10.0,
    )

    # Build grain boundary with vacuum
    builder = GrainBoundaryWithVacuumBuilder()
    grain_boundary_slab = builder.get_material(config)

    assert isinstance(grain_boundary_slab, Material)
    assert len(grain_boundary_slab.basis.elements.values) > 0
