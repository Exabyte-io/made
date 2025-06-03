from mat3ra.made.material import Material
from mat3ra.made.tools.build.slab import SlabBuilderParameters, SlabConfiguration, create_slab
from mat3ra.made.tools.build.slab.configuration import (
    AtomicLayersUniqueRepeated,
    Termination,
    VacuumConfiguration,
    get_terminations,
    choose_termination,
)
from unit.fixtures.slab import SI_SLAB_001, SI_SLAB_001_CONFIGURATION, SI_SLAB_DEFAULT_PARAMETERS

from .utils import assert_two_entities_deep_almost_equal

MILLER_INDICES = SI_SLAB_001_CONFIGURATION["miller_indices"]
USE_CONVENTIONAL_CELL = SI_SLAB_001_CONFIGURATION["use_conventional_cell"]
NUMBER_OF_LAYERS = SI_SLAB_001_CONFIGURATION["number_of_layers"]
VACUUM = SI_SLAB_001_CONFIGURATION["vacuum"]
XY_SUPERCELL_MATRIX = SI_SLAB_001_CONFIGURATION["xy_supercell_matrix"]

material = Material.create_default()


def test_build_slab_human_usage():
    crystal = material
    hkl = MILLER_INDICES
    terminations = get_terminations(crystal=crystal, miller_indices=hkl)
    termination = choose_termination(terminations, "Si")
    configuration = SlabConfiguration.from_parameters(
        crystal=crystal,
        miller_indices=hkl,
        use_conventional_cell=USE_CONVENTIONAL_CELL,
        termination=termination,
        number_of_layers=NUMBER_OF_LAYERS,
        vacuum=VACUUM,
        xy_supercell_matrix=XY_SUPERCELL_MATRIX,
    )

    slab = create_slab(configuration=configuration, build_parameters=SlabBuilderParameters())
    assert_two_entities_deep_almost_equal(slab, SI_SLAB_001)
