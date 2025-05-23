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
from unit.fixtures.slab import SI_SLAB_100, SI_SLAB_DEFAULT_PARAMETERS

from .utils import assert_two_entities_deep_almost_equal

material = Material.create_default()


def test_get_terminations():

    crystal_lattice_planes = CrystalLatticePlanes(crystal=material, miller_indices=[0, 0, 1])
    terminations = crystal_lattice_planes.get_terminations()
    print("Available terminations:", [str(t) for t in terminations])
    assert len(terminations) > 0, "No terminations found for the given crystal and miller indices."


def test_build_slab():
    crystal_lattice_planes = CrystalLatticePlanes(crystal=material, miller_indices=[0, 0, 1])
    terminations = crystal_lattice_planes.get_terminations()
    print("Available terminations:", [str(t) for t in terminations])

    # Use the available termination for both the atomic layers and the reference
    atomic_layers = AtomicLayersUniqueRepeated(
        crystal=material,
        miller_indices=[0, 0, 1],
        number_of_repetitions=2,
        termination_top=terminations[0],
        use_conventional_cell=True,
    )

    vacuum = VacuumConfiguration(direction=AxisEnum.z, size=5.0, is_orthogonal=True)
    slab_config = SlabConfiguration(
        supercell_xy=[[1, 0], [0, 1]], stack_components=[atomic_layers, vacuum], direction=AxisEnum.z
    )
    params = SlabBuilderParameters(min_vacuum_size=5.0, reorient_lattice=True, symmetrize=True)
    slab = create_slab(slab_config, build_parameters=params)
    # Use SI_SLAB_100 as reference, but update if needed to match the new termination
    assert str(terminations[0]) == "Si_P4/mmm_2", "Test expects Si_P4/mmm_2 termination."
    assert_two_entities_deep_almost_equal(slab, SI_SLAB_100)


def test_build_slab_with_default_parameters():
    # Create atomic layers component with default parameters
    atomic_layers = AtomicLayersUniqueRepeated(
        crystal=material,
        miller_indices=[0, 0, 1],
        number_of_repetitions=1,
        termination_top=Termination.from_string("Si_P4/mmm_1"),
        use_conventional_cell=True,
    )

    # Create vacuum component with default parameters
    vacuum = VacuumConfiguration(direction=AxisEnum.z)

    # Create slab configuration with required fields
    slab_config = SlabConfiguration(
        supercell_xy=[[1, 0], [0, 1]], stack_components=[atomic_layers, vacuum], direction=AxisEnum.z
    )
    slab = create_slab(slab_config)
    assert_two_entities_deep_almost_equal(slab, SI_SLAB_DEFAULT_PARAMETERS)
