from mat3ra.made.material import Material
from mat3ra.made.tools.build.slab import (
    SlabBuilderParameters,
    SlabConfiguration,
    create_slab,
    get_terminations,
)
from mat3ra.made.tools.build.slab.configuration import (
    CrystalLatticePlanes,
    AtomicLayersUniqueRepeated,
    VacuumConfiguration,
    Termination,
)
from mat3ra.esse.models.materials_category.pristine_structures.two_dimensional.slab import AxisEnum
from unit.fixtures.slab import SI_SLAB_001, SI_SLAB_DEFAULT_PARAMETERS, SI_SLAB_001_2_ATOMS, SI_SLAB_001_CONFIGURATION

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
    expected_termination = Termination.from_string("Si_R-3m_1")
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
    vacuum = VacuumConfiguration(direction=AxisEnum.z, size=VACUUM)
    slab_config = SlabConfiguration(
        xy_supercell_matrix=XY_SUPERCELL_MATRIX, stack_components=[atomic_layers, vacuum], direction=AxisEnum.z
    )
    params = SlabBuilderParameters(min_vacuum_size=1, reorient_lattice=True, symmetrize=True, make_primitive=True)

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
    # Create atomic layers component with default parameters
    atomic_layers = AtomicLayersUniqueRepeated(
        crystal=material,
        miller_indices=MILLER_INDICES,
        number_of_repetitions=1,
        termination_top=Termination.from_string("Si_P4/mmm_1"),
        use_conventional_cell=True,
    )

    # Create vacuum component with default parameters
    vacuum = VacuumConfiguration(direction=AxisEnum.z)

    # Create slab configuration with required fields
    slab_config = SlabConfiguration(
        xy_supercell_matrix=[[1, 0], [0, 1]], stack_components=[atomic_layers, vacuum], direction=AxisEnum.z
    )
    slab = create_slab(slab_config)
    assert_two_entities_deep_almost_equal(slab, SI_SLAB_DEFAULT_PARAMETERS)
