import sys

import pytest
from mat3ra.made.lattice import COORDINATE_TOLERANCE
from mat3ra.made.material import Material
from mat3ra.made.tools.build.defect import EquidistantAdatomSlabDefectBuilder, create_slab_defect
from mat3ra.made.tools.build.defect.builders import CrystalSiteAdatomSlabDefectBuilder
from mat3ra.made.tools.build.defect.configuration import AdatomSlabPointDefectConfiguration
from mat3ra.utils import assertion as assertion_utils

from mat3ra.made.tools.build.defect.enums import AdatomPlacementMethodEnum
from mat3ra.made.tools.build.defect.slab.helpers import create_adatom_defect
from unit.fixtures.slab import SI_CONVENTIONAL_SLAB_001


@pytest.mark.parametrize(
    "crystal_config, position_on_surface, distance_z, chemical_element, adatom_placement_method, expected_last_coord",
    [
        (
            SI_CONVENTIONAL_SLAB_001,
            [0.5, 0.5],
            2.0,
            "Si",
            AdatomPlacementMethodEnum.NEW_CRYSTAL_SITE,
            [0.5, 0.5, 0.57482],
        ),
        (
            SI_CONVENTIONAL_SLAB_001,
            [0.55, 0.51],
            2.0,
            "Si",
            AdatomPlacementMethodEnum.EQUIDISTANT,
            [0.5, 0.5, 0.72598],
        ),
    ],
)
def test_create_adatom(
    crystal_config,
    position_on_surface,
    distance_z,
    chemical_element,
    adatom_placement_method,
    expected_last_coord,
):
    slab = Material.create(crystal_config)
    defect = create_adatom_defect(
        slab, position_on_surface, distance_z, adatom_placement_method, chemical_element, distance_z
    )
    assertion_utils.assert_deep_almost_equal(expected_last_coord, defect.basis.coordinates.values[-1])


@pytest.mark.parametrize(
    "slab_material_config, position_on_surface, distance_z, expected_center",
    [(SI_CONVENTIONAL_SLAB_001, [0.5432, 0.5123], 2.5, [0.5, 0.5, 0.75735])],
)
def test_get_equidistant_position(slab_material_config, position_on_surface, distance_z, expected_center):
    builder = EquidistantAdatomSlabDefectBuilder()

    slab_material = Material.create(slab_material_config)

    equidistant_position = builder.get_equidistant_position(
        material=slab_material, position_on_surface=position_on_surface, distance_z=distance_z
    )

    assertion_utils.assert_deep_almost_equal(equidistant_position, expected_center)
