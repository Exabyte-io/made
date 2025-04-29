from mat3ra.made.material import Material
from mat3ra.made.tools.build.slab import (
    PymatgenSlabGeneratorParameters,
    SlabConfiguration,
    create_slab,
    get_terminations,
)
from unit.fixtures.slab import SI_SLAB_100

from .utils import assert_two_entities_deep_almost_equal

material = Material.create_default()


def test_build_slab():
    slab_config = SlabConfiguration(
        bulk=material,
        miller_indices=(0, 0, 1),
        thickness=2,
        vacuum=5.0,
        xy_supercell_matrix=[[1, 0], [0, 1]],
        use_orthogonal_z=True,
        use_conventional_cell=True,
        make_primitive=True,
    )
    params = PymatgenSlabGeneratorParameters(
        min_vacuum_size=1, in_unit_planes=True, reorient_lattice=True, symmetrize=True
    )
    terminations = get_terminations(slab_config, build_parameters=params)
    termination = terminations[0]
    slab = create_slab(slab_config, termination, params)
    assert_two_entities_deep_almost_equal(slab, SI_SLAB_100)
