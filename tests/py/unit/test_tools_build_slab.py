from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.lattice_planes import CrystalLatticePlanesMaterialAnalyzer
from mat3ra.made.tools.build.slab import (
    SlabConfiguration,
    create_slab,
    SlabBuilder,
    select_slab_termination,
    get_slab_terminations,
)
from mat3ra.made.tools.build.slab.builders import (
    AtomicLayersUniqueRepeatedBuilder,
    SlabBuilderParameters,
)
from mat3ra.made.tools.build.slab.configuration import (
    AtomicLayersUniqueRepeatedConfiguration,
    VacuumConfiguration,
)
from unit.fixtures.slab import SI_SLAB_001, SI_SLAB_001_CONFIGURATION
from .utils import assert_two_entities_deep_almost_equal

MILLER_INDICES = SI_SLAB_001_CONFIGURATION["miller_indices"]
USE_CONVENTIONAL_CELL = SI_SLAB_001_CONFIGURATION["use_conventional_cell"]
NUMBER_OF_LAYERS = SI_SLAB_001_CONFIGURATION["number_of_layers"]
VACUUM = SI_SLAB_001_CONFIGURATION["vacuum"]
XY_SUPERCELL_MATRIX = SI_SLAB_001_CONFIGURATION["xy_supercell_matrix"]

material = Material.create_default()


def test_build_slab():
    crystal_lattice_planes_analyzer = CrystalLatticePlanesMaterialAnalyzer(
        material=material, miller_indices=MILLER_INDICES
    )
    terminations = crystal_lattice_planes_analyzer.terminations
    termination = select_slab_termination(terminations, "Si")

    atomic_layers_repeated_configuration = AtomicLayersUniqueRepeatedConfiguration(
        crystal=material,
        miller_indices=MILLER_INDICES,
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

    assert_two_entities_deep_almost_equal(slab, SI_SLAB_001)


def test_create_slab():
    crystal = material
    hkl = MILLER_INDICES
    terminations = get_slab_terminations(crystal=crystal, miller_indices=hkl)
    termination = select_slab_termination(terminations, "Si")
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
