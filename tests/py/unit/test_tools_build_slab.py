from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from mat3ra.made.material import Material
from mat3ra.made.tools.build.slab import (
    SlabBuilderParameters,
    SlabConfiguration,
    create_slab,
)
from mat3ra.made.tools.build.slab.configuration import (
    CrystalLatticePlanes,
    AtomicLayersUniqueRepeated,
    VacuumConfiguration,
    Termination,
)
from unit.fixtures.slab import SI_SLAB_001, SI_SLAB_DEFAULT_PARAMETERS, SI_SLAB_001_CONFIGURATION
from .utils import assert_two_entities_deep_almost_equal

material = Material.create_default()
MILLER_INDICES = SI_SLAB_001_CONFIGURATION["miller_indices"]
USE_CONVENTIONAL_CELL = SI_SLAB_001_CONFIGURATION["use_conventional_cell"]
NUMBER_OF_LAYERS = SI_SLAB_001_CONFIGURATION["number_of_layers"]
VACUUM = SI_SLAB_001_CONFIGURATION["vacuum"]
XY_SUPERCELL_MATRIX = SI_SLAB_001_CONFIGURATION["xy_supercell_matrix"]


def test_get_terminations():
    crystal_lattice_planes = CrystalLatticePlanes(crystal=material, miller_indices=MILLER_INDICES)
    terminations = crystal_lattice_planes.get_terminations()
    print("Available terminations:", [str(t) for t in terminations])
    expected_termination = Termination.from_string("Si_Fm-3m_1")
    assert_two_entities_deep_almost_equal(terminations[0], expected_termination)


def test_build_slab():
    crystal_lattice_planes = CrystalLatticePlanes(
        crystal=material, miller_indices=MILLER_INDICES, use_conventional_cell=USE_CONVENTIONAL_CELL
    )
    terminations = crystal_lattice_planes.get_terminations()
    print("Available terminations:", [str(t) for t in terminations])

    # Use the available termination for both the atomic layers and the reference
    atomic_layers = AtomicLayersUniqueRepeated(
        crystal=material,
        miller_indices=MILLER_INDICES,
        use_conventional_cell=USE_CONVENTIONAL_CELL,
        number_of_repetitions=NUMBER_OF_LAYERS,
        termination_top=terminations[0],
    )
    vacuum = VacuumConfiguration(direction="z", size=VACUUM)
    slab_config = SlabConfiguration(
        xy_supercell_matrix=XY_SUPERCELL_MATRIX, stack_components=[atomic_layers, vacuum], direction="z"
    )
    params = SlabBuilderParameters(min_vacuum_size=0, reorient_lattice=True, symmetrize=True, make_primitive=False)

    slab = create_slab(slab_config, build_parameters=params)

    assert_two_entities_deep_almost_equal(slab, SI_SLAB_001)


def test_build_slab_from_parameters():
    slab_config = SlabConfiguration.from_parameters(
        bulk=material,
        miller_indices=MILLER_INDICES,
        number_of_layers=NUMBER_OF_LAYERS,
        vacuum=VACUUM,
        use_conventional_cell=USE_CONVENTIONAL_CELL,
    )
    slab = create_slab(slab_config)
    assert_two_entities_deep_almost_equal(slab, SI_SLAB_001)


def test_build_slab_with_default_parameters():

    # Create slab configuration with minimal set up
    slab_config = SlabConfiguration.from_parameters(
        bulk=material,
    )
    slab = create_slab(slab_config)
    assert_two_entities_deep_almost_equal(slab, SI_SLAB_DEFAULT_PARAMETERS)
