import copy
from typing import Final, Tuple

import numpy as np
import pytest
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.lattice_planes import CrystalLatticePlanesMaterialAnalyzer
from mat3ra.made.tools.build import MaterialWithBuildMetadata
from mat3ra.made.tools.build.slab.builders import AtomicLayersUniqueRepeatedBuilder, SlabBuilder, SlabBuilderParameters
from mat3ra.made.tools.build.slab.configurations import AtomicLayersUniqueRepeatedConfiguration, SlabConfiguration
from mat3ra.made.tools.build.slab.helpers import create_slab, get_slab_terminations
from mat3ra.made.tools.build.slab.termination_utils import select_slab_termination
from mat3ra.made.tools.build.vacuum.configuration import VacuumConfiguration
from mat3ra.made.utils import AXIS_TO_INDEX_MAP, adjust_material_cell_to_set_gap_along_direction
from mat3ra.utils import assertion
from unit.fixtures.bulk import BULK_Si_CONVENTIONAL, BULK_Si_PRIMITIVE
from unit.fixtures.slab import (
    SI_CONVENTIONAL_SLAB_001,
    SI_PRIMITIVE_SLAB_001,
    SLAB_SI_CONVENTIONAL_001_NO_VACUUM,
    SLAB_SrTiO3_011_TERMINATION_O2,
    SLAB_SrTiO3_011_TERMINATION_SrTiO,
)

from .fixtures.generated.fixtures import BULK_SrTiO3
from .fixtures.monolayer import GRAPHENE
from .utils import assert_two_entities_deep_almost_equal

SI_CONVENTIONAL_SLAB_001_NO_BUILD_METADATA = copy.deepcopy(SI_CONVENTIONAL_SLAB_001)
SI_CONVENTIONAL_SLAB_001_NO_BUILD_METADATA["metadata"].pop("build", None)

PARAMS_SI_PRIMITIVE_SLAB_001: Final = (
    BULK_Si_PRIMITIVE,
    (0, 0, 1),
    "Si",
    2,
    5.0,
    [[1, 0], [0, 1]],
)

PARAMS_BUILD_SLAB_CONVENTIONAL: Final = (
    BULK_Si_PRIMITIVE,
    (0, 0, 1),
    "Si",
    2,
    5.0,
    [[1, 0], [0, 1]],
)

PARAMS_BUILD_SLAB_CONVENTIONAL_NO_VACUUM: Final = (
    BULK_Si_PRIMITIVE,
    (0, 0, 1),
    "Si",
    2,
    0.0,
    [[1, 0], [0, 1]],
)

PARAMS_BUILD_SLAB_CONVENTIONAL_SrTiO_SrTiO: Final = (
    BULK_SrTiO3,
    (0, 1, 1),
    "SrTiO",
    2,
    5.0,
    [[1, 0], [0, 1]],
)

PARAMS_BUILD_SLAB_CONVENTIONAL_SrTiO_O2: Final = (
    BULK_SrTiO3,
    (0, 1, 1),
    "O2",
    2,
    5.0,
    [[1, 0], [0, 1]],
)

PARAMS_CREATE_SLAB: Final = (
    BULK_Si_CONVENTIONAL,
    (0, 0, 1),
    "Si",
    2,
    5,
    [[1, 0], [0, 1]],
    True,
)


def get_slab_with_builder(
    material: Material,
    miller_indices: Tuple[int, int, int],
    termination_formula: str,
    number_of_layers: int,
    vacuum: float,
    xy_supercell_matrix: list,
) -> Material:
    crystal_lattice_planes_analyzer = CrystalLatticePlanesMaterialAnalyzer(
        material=material, miller_indices=miller_indices
    )
    terminations = crystal_lattice_planes_analyzer.terminations
    termination = select_slab_termination(terminations, termination_formula)

    atomic_layers_repeated_configuration = AtomicLayersUniqueRepeatedConfiguration(
        crystal=material,
        miller_indices=miller_indices,
        termination_top=termination,
        number_of_repetitions=number_of_layers,
    )
    atomic_layers_repeated_orthogonal_c = AtomicLayersUniqueRepeatedBuilder().get_material(
        atomic_layers_repeated_configuration
    )
    vacuum_configuration = VacuumConfiguration(
        size=vacuum, crystal=atomic_layers_repeated_orthogonal_c, direction=AxisEnum.z
    )
    build_params = SlabBuilderParameters(use_orthogonal_c=True, xy_supercell_matrix=xy_supercell_matrix)
    slab_configuration = SlabConfiguration(
        stack_components=[atomic_layers_repeated_configuration, vacuum_configuration],
        direction=AxisEnum.z,
    )
    builder = SlabBuilder(build_parameters=build_params)
    slab = builder.get_material(slab_configuration)

    return slab


@pytest.mark.parametrize(
    "material_config, miller_indices, termination_formula, number_of_layers,"
    + " vacuum, xy_supercell_matrix, expected_slab_config",
    [
        (*PARAMS_SI_PRIMITIVE_SLAB_001, SI_PRIMITIVE_SLAB_001),
    ],
)
def test_build_slab_primitive(
    material_config,
    miller_indices,
    termination_formula,
    number_of_layers,
    vacuum,
    xy_supercell_matrix,
    expected_slab_config,
):
    material = MaterialWithBuildMetadata.create(material_config)
    slab = get_slab_with_builder(
        material, miller_indices, termination_formula, number_of_layers, vacuum, xy_supercell_matrix
    )
    slab.metadata.build = []  # Remove build metadata for comparison
    expected_slab_config.get("metadata", {}).pop("build", None)  # Remove build metadata for comparison
    assert_two_entities_deep_almost_equal(slab, expected_slab_config)


@pytest.mark.parametrize(
    "material_config, miller_indices, termination_formula, number_of_layers,"
    + " vacuum, xy_supercell_matrix, expected_slab_config",
    [
        (
            *PARAMS_BUILD_SLAB_CONVENTIONAL,
            SI_CONVENTIONAL_SLAB_001_NO_BUILD_METADATA,
        ),
        (
            *PARAMS_BUILD_SLAB_CONVENTIONAL_NO_VACUUM,
            SLAB_SI_CONVENTIONAL_001_NO_VACUUM,
        ),
    ],
)
def test_build_slab_conventional(
    material_config,
    miller_indices,
    termination_formula,
    number_of_layers,
    vacuum,
    xy_supercell_matrix,
    expected_slab_config,
):
    material = Material.create(material_config)
    crystal_lattice_planes_analyzer = CrystalLatticePlanesMaterialAnalyzer(
        material=material, miller_indices=miller_indices
    )
    conventional_material = crystal_lattice_planes_analyzer.material_with_conventional_lattice
    slab = get_slab_with_builder(
        conventional_material,
        miller_indices,
        termination_formula,
        number_of_layers,
        vacuum,
        xy_supercell_matrix,
    )
    slab.metadata.build = []
    assert_two_entities_deep_almost_equal(slab, expected_slab_config)


@pytest.mark.parametrize(
    "material_config, miller_indices, termination_formula, number_of_layers,"
    + " vacuum, xy_supercell_matrix, expected_slab_config",
    [
        (
            *PARAMS_BUILD_SLAB_CONVENTIONAL_SrTiO_SrTiO,
            SLAB_SrTiO3_011_TERMINATION_SrTiO,
        ),
        (
            *PARAMS_BUILD_SLAB_CONVENTIONAL_SrTiO_O2,
            SLAB_SrTiO3_011_TERMINATION_O2,
        ),
    ],
)
def test_build_slab_conventional_with_multiple_terminations(
    material_config,
    miller_indices,
    termination_formula,
    number_of_layers,
    vacuum,
    xy_supercell_matrix,
    expected_slab_config,
):
    crystal_lattice_planes_analyzer = CrystalLatticePlanesMaterialAnalyzer(
        material=material_config, miller_indices=miller_indices
    )
    conventional_material = crystal_lattice_planes_analyzer.material_with_conventional_lattice
    slab = get_slab_with_builder(
        conventional_material,
        miller_indices,
        termination_formula,
        number_of_layers,
        vacuum,
        xy_supercell_matrix,
    )

    slab.metadata.build = []  # Remove build metadata for comparison
    expected_slab_config.get("metadata", {}).pop("build", None)
    assert_two_entities_deep_almost_equal(slab, expected_slab_config)


@pytest.mark.parametrize(
    "material_config, miller_indices, termination_formula, number_of_layers,"
    + " vacuum, xy_supercell, use_conventional_cell, expected_slab_config",
    [
        (
            *PARAMS_CREATE_SLAB,
            SI_CONVENTIONAL_SLAB_001_NO_BUILD_METADATA,
        ),
    ],
)
def test_create_slab(
    material_config,
    miller_indices,
    termination_formula,
    number_of_layers,
    vacuum,
    xy_supercell,
    use_conventional_cell,
    expected_slab_config,
):
    crystal = Material.create(material_config)
    terminations = get_slab_terminations(material=crystal, miller_indices=miller_indices)
    termination = select_slab_termination(terminations, termination_formula)
    slab = create_slab(
        crystal=crystal,
        miller_indices=miller_indices,
        use_conventional_cell=use_conventional_cell,
        termination=termination,
        number_of_layers=number_of_layers,
        vacuum=vacuum,
        xy_supercell_matrix=xy_supercell,
    )
    slab.metadata.build = []  # Remove build metadata for comparison
    assert_two_entities_deep_almost_equal(slab, expected_slab_config)


@pytest.mark.parametrize(
    "material_config, direction, gap, expected_length",
    [
        (SLAB_SI_CONVENTIONAL_001_NO_VACUUM, AxisEnum.z, 5.0, 14.5703),  # Adjusted length for z direction
        (GRAPHENE, AxisEnum.y, 5.0, 6.6448),  # Adjusted length for y direction (lattice vector b)
    ],
)
def test_adjust_lattice_for_gap(material_config, direction, gap, expected_length):
    material = Material.create(material_config)
    adjusted_material = adjust_material_cell_to_set_gap_along_direction(material, gap, direction)

    axis_index = AXIS_TO_INDEX_MAP[direction.value]
    actual_length = np.linalg.norm(adjusted_material.lattice.vector_arrays[axis_index])
    assertion.assert_almost_equal_numbers(actual_length.tolist(), expected_length)
