from typing import Tuple

# import pytest
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.lattice_planes import CrystalLatticePlanesMaterialAnalyzer
from mat3ra.made.tools.build.slab.builders import AtomicLayersUniqueRepeatedBuilder, SlabBuilder, SlabBuilderParameters
from mat3ra.made.tools.build.slab.configuration import (
    AtomicLayersUniqueRepeatedConfiguration,
    SlabConfiguration,
    VacuumConfiguration,
)
from mat3ra.made.tools.build.slab.helpers import create_slab, get_slab_terminations, select_slab_termination
from unit.fixtures.bulk import SI_PRIMITIVE_CELL_MATERIAL
from unit.fixtures.slab import (
    SI_CONVENTIONAL_SLAB_001,
    SI_PRIMITIVE_SLAB_001,
    SI_SLAB_001_CONFIGURATION_FROM_CONVENTIONAL,
)

from .utils import assert_two_entities_deep_almost_equal

MILLER_INDICES = SI_SLAB_001_CONFIGURATION_FROM_CONVENTIONAL["miller_indices"]
USE_CONVENTIONAL_CELL = SI_SLAB_001_CONFIGURATION_FROM_CONVENTIONAL["use_conventional_cell"]
NUMBER_OF_LAYERS = SI_SLAB_001_CONFIGURATION_FROM_CONVENTIONAL["number_of_layers"]
VACUUM = SI_SLAB_001_CONFIGURATION_FROM_CONVENTIONAL["vacuum"]
XY_SUPERCELL_MATRIX = SI_SLAB_001_CONFIGURATION_FROM_CONVENTIONAL["xy_supercell_matrix"]

SI_CONVENTIONAL_SLAB_001_NO_BUILD_METADATA = SI_CONVENTIONAL_SLAB_001.copy()
SI_CONVENTIONAL_SLAB_001_NO_BUILD_METADATA["metadata"].pop("build", None)


def get_slab_with_builder(
    material: Material, miller_indices: Tuple[int, int, int], termination_formula: str
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
        number_of_repetitions=NUMBER_OF_LAYERS,
    )
    atomic_layers_repeated_orthogonal_c = AtomicLayersUniqueRepeatedBuilder().get_material(
        atomic_layers_repeated_configuration
    )
    vacuum_configuration = VacuumConfiguration(
        size=VACUUM, crystal=atomic_layers_repeated_orthogonal_c, direction=AxisEnum.z
    )
    build_params = SlabBuilderParameters(use_orthogonal_c=True, xy_supercell_matrix=XY_SUPERCELL_MATRIX)
    slab_configuration = SlabConfiguration(
        stack_components=[atomic_layers_repeated_configuration, vacuum_configuration],
        direction=AxisEnum.z,
    )
    builder = SlabBuilder(build_parameters=build_params)
    slab = builder.get_material(slab_configuration)

    return slab


def test_build_slab_primitive():
    slab = get_slab_with_builder(SI_PRIMITIVE_CELL_MATERIAL, MILLER_INDICES, "Si")
    slab.metadata.pop("build")  # Remove build metadata for comparison
    SI_PRIMITIVE_SLAB_001["metadata"].pop("build")  # Remove build metadata for comparison
    assert_two_entities_deep_almost_equal(slab, SI_PRIMITIVE_SLAB_001)


def test_build_slab_conventional():
    crystal_lattice_planes_analyzer = CrystalLatticePlanesMaterialAnalyzer(
        material=SI_PRIMITIVE_CELL_MATERIAL, miller_indices=MILLER_INDICES
    )
    conventional_material = crystal_lattice_planes_analyzer.material_with_conventional_lattice
    slab = get_slab_with_builder(conventional_material, MILLER_INDICES, "Si")
    slab.metadata.pop("build")
    assert_two_entities_deep_almost_equal(slab, SI_CONVENTIONAL_SLAB_001_NO_BUILD_METADATA)


def test_build_slab_conventional_with_multiple_terminations():
    from mat3ra.standata.materials import Materials

    SrTiO_BULK = Material.create(Materials().get_by_name_first_match("SrTiO3"))
    SrTiO_MILLER_INDICES = (0, 1, 1)
    crystal_lattice_planes_analyzer = CrystalLatticePlanesMaterialAnalyzer(
        material=SrTiO_BULK, miller_indices=SrTiO_MILLER_INDICES
    )
    conventional_material = crystal_lattice_planes_analyzer.material_with_conventional_lattice
    terminations = get_slab_terminations(material=conventional_material, miller_indices=SrTiO_MILLER_INDICES)
    termination_1 = select_slab_termination(terminations, "SrTiO")
    slab_1 = get_slab_with_builder(conventional_material, SrTiO_MILLER_INDICES, termination_1.formula)
    termination_2 = select_slab_termination(terminations, "O2")
    slab_2 = get_slab_with_builder(conventional_material, SrTiO_MILLER_INDICES, termination_2.formula)

    slab_1.metadata.pop("build")  # Remove build metadata for comparison
    slab_2.metadata.pop("build")  # Remove build metadata for comparison
    assert_two_entities_deep_almost_equal(slab_1, SI_CONVENTIONAL_SLAB_001_NO_BUILD_METADATA)


def test_create_slab():
    crystal = SI_PRIMITIVE_CELL_MATERIAL
    terminations = get_slab_terminations(material=crystal, miller_indices=MILLER_INDICES)
    termination = select_slab_termination(terminations, "Si")
    slab = create_slab(
        crystal=crystal,
        miller_indices=MILLER_INDICES,
        use_conventional_cell=USE_CONVENTIONAL_CELL,
        termination=termination,
        number_of_layers=NUMBER_OF_LAYERS,
        vacuum=VACUUM,
        xy_supercell_matrix=XY_SUPERCELL_MATRIX,
    )
    slab.metadata.pop("build")  # Remove build metadata for comparison
    assert_two_entities_deep_almost_equal(slab, SI_CONVENTIONAL_SLAB_001_NO_BUILD_METADATA)
