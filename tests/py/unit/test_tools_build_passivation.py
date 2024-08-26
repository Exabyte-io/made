from mat3ra.made.material import Material

from mat3ra.made.tools.modify import (
    passivate_surface,
)
from mat3ra.utils import assertion as assertion_utils

from .fixtures import SI_SLAB, SI_SLAB_PASSIVATED


def test_passivate_surface():
    passivated_material = passivate_surface(
        material=Material(SI_SLAB), passivant="H", default_bond_length=1.48, surface="both"
    )
    assertion_utils.assert_deep_almost_equal(SI_SLAB_PASSIVATED, passivated_material.to_json())
