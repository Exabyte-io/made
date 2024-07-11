from mat3ra.made.material import Material
from mat3ra.made.tools.build.slab import SlabConfiguration, create_slab, get_terminations
from mat3ra.utils import assertion as assertion_utils

from .fixtures import SI_SLAB

material = Material.create(Material.default_config)


def test_build_slab():
    slab_config = SlabConfiguration(
        bulk=material,
        miller_indices=(0, 0, 1),
        thickness=1,
        vacuum=1,
        xy_supercell_matrix=[[1, 0], [0, 1]],
        use_orthogonal_z=True,
        make_primitive=True,
    )
    termination = get_terminations(slab_config)[0]
    slab = create_slab(slab_config, termination)
    assertion_utils.assert_deep_almost_equal(SI_SLAB, slab.to_json())
