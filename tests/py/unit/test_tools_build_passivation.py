from mat3ra.made.material import Material
from mat3ra.made.tools.build.passivation.builders import SurfacePassivationBuilder, SurfacePassivationBuilderParameters
from mat3ra.made.tools.build.passivation.configuration import PassivationConfiguration
from mat3ra.utils import assertion as assertion_utils

from .fixtures import SI_SLAB, SI_SLAB_PASSIVATED


def test_passivate_surface():
    config = PassivationConfiguration(slab=Material(SI_SLAB), passivant="H", bond_length=1.48, surface="both")
    builder = SurfacePassivationBuilder(
        build_parameters=SurfacePassivationBuilderParameters(shadowing_radius=2.5, depth=2.0)
    )
    passivated_material = builder.get_material(config)
    assertion_utils.assert_deep_almost_equal(SI_SLAB_PASSIVATED, passivated_material.to_json())
