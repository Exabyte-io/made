from mat3ra.made.material import Material
from mat3ra.made.tools.build.slab import (
    PymatgenSlabGeneratorParameters,
    SlabConfiguration,
    create_slab,
    get_terminations,
)
from mat3ra.utils import assertion as assertion_utils

from .fixtures import SI_SLAB_100_MAVRL

material = Material.create(Material.default_config)


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
    termination = get_terminations(slab_config, build_parameters=params)[0]
    slab = create_slab(slab_config, termination, params)
    print(slab.to_json())
    assertion_utils.assert_deep_almost_equal(SI_SLAB_100_MAVRL, slab.to_json())
