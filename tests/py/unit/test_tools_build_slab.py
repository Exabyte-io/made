from mat3ra.code.vector import Vector3D
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from mat3ra.made.material import Material
from mat3ra.made.tools.build.slab import SlabConfiguration, create_slab, SlabBuilder
from mat3ra.made.tools.build.slab.configuration import (
    AtomicLayersUniqueRepeatedConfiguration,
    VacuumConfiguration,
    get_terminations,
    choose_termination,
    MillerSupercell,
    CrystalLatticePlanesConfiguration,
    AtomicLayersUnique,
    CrystalLatticePlanesBuilder,
    AtomicLayersUniqueRepeatedBuilder,
    ConventionalCellConfiguration,
    ConventionalCellBuilder,
)
from mat3ra.made.tools.operations.core.unary import translate
from unit.fixtures.slab import SI_SLAB_001, SI_SLAB_001_CONFIGURATION
from .utils import assert_two_entities_deep_almost_equal

MILLER_INDICES = SI_SLAB_001_CONFIGURATION["miller_indices"]
USE_CONVENTIONAL_CELL = SI_SLAB_001_CONFIGURATION["use_conventional_cell"]
NUMBER_OF_LAYERS = SI_SLAB_001_CONFIGURATION["number_of_layers"]
VACUUM = SI_SLAB_001_CONFIGURATION["vacuum"]
XY_SUPERCELL_MATRIX = SI_SLAB_001_CONFIGURATION["xy_supercell_matrix"]

material = Material.create_default()


def test_build_slab():
    conventional_cell_config = ConventionalCellConfiguration(crystal=material)
    conventional_cell = (
        ConventionalCellBuilder().get_material(conventional_cell_config) if USE_CONVENTIONAL_CELL else material
    )
    crystal_lattice_planes = CrystalLatticePlanesConfiguration(crystal=conventional_cell, miller_indices=MILLER_INDICES)
    terminations = crystal_lattice_planes.terminations
    termination = choose_termination(terminations, "Si")

    atomic_layers = AtomicLayersUnique(
        crystal=crystal_lattice_planes.crystal,
        miller_indices=MILLER_INDICES,
    )

    atomic_layers_repeated_config = AtomicLayersUniqueRepeatedConfiguration(
        crystal=crystal_lattice_planes.crystal,
        miller_indices=MILLER_INDICES,
        termination_top=termination,
        number_of_repetitions=NUMBER_OF_LAYERS,
    )
    atomic_layers_repeated_orthogonal_c = AtomicLayersUniqueRepeatedBuilder().get_material(
        atomic_layers_repeated_config
    )

    translation_vector: Vector3D = atomic_layers.get_translation_vector(termination)
    atomic_layers_repeated_terminated = translate(atomic_layers_repeated_orthogonal_c, translation_vector)

    vacuum_configuration = VacuumConfiguration(
        size=VACUUM, crystal=atomic_layers_repeated_terminated, direction=AxisEnum.z
    )

    slab_configuration = SlabConfiguration(
        stack_components=[atomic_layers_repeated_terminated, vacuum_configuration],
        direction=AxisEnum.z,
        xy_supercell_matrix=XY_SUPERCELL_MATRIX,
    )

    builder = SlabBuilder()
    slab = builder.get_material(slab_configuration)

    assert_two_entities_deep_almost_equal(slab, SI_SLAB_001)


def test_build_slab_human_usage():
    crystal = material
    hkl = MILLER_INDICES
    crystal_planes_configuration = CrystalLatticePlanesConfiguration(
        crystal=crystal,
        miller_indices=hkl,
        use_conventional_cell=USE_CONVENTIONAL_CELL,
    )

    terminations = get_terminations(crystal=crystal, miller_indices=hkl)
    crystal_planes = CrystalLatticePlanesBuilder().get_material(crystal_planes_configuration)
    # visualize(crystal_planes)
    termination = choose_termination(terminations, "Si")
    slab = create_slab(
        crystal=crystal,
        miller_indices=hkl,
        use_conventional_cell=USE_CONVENTIONAL_CELL,
        termination=termination,
        number_of_layers=NUMBER_OF_LAYERS,
        vacuum=VACUUM,
        xy_supercell_matrix=XY_SUPERCELL_MATRIX,
    )
    assert_two_entities_deep_almost_equal(slab, SI_SLAB_001)
