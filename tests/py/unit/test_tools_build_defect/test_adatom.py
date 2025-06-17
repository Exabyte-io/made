import sys

import pytest
from mat3ra.made.lattice import COORDINATE_TOLERANCE
from mat3ra.made.material import Material
from mat3ra.made.tools.build.defect import EquidistantAdatomSlabDefectBuilder, create_slab_defect
from mat3ra.made.tools.build.defect.builders import CrystalSiteAdatomSlabDefectBuilder
from mat3ra.made.tools.build.defect.configuration import AdatomSlabPointDefectConfiguration
from mat3ra.utils import assertion as assertion_utils
from unit.fixtures.slab import SI_CONVENTIONAL_SLAB_001


@pytest.mark.skip(reason="we'll fix before epic-7623 is merged")
@pytest.mark.parametrize(
    "crystal_config, position_on_surface, distance_z, chemical_element, expected_last_element, expected_last_coord",
    [
        (
            SI_CONVENTIONAL_SLAB_001,
            [0.5, 0.5],
            2,
            "Si",
            "Si",
            [0.5, 0.5, 0.7641022175421057],
        )
    ],
)
def test_create_adatom(
    crystal_config, position_on_surface, distance_z, chemical_element, expected_last_element, expected_last_coord
):
    # Adatom of Si at 0.5, 0.5 position
    crystal = Material.create(crystal_config)
    configuration = AdatomSlabPointDefectConfiguration(
        crystal=crystal,
        position_on_surface=position_on_surface,
        distance_z=distance_z,
        chemical_element=chemical_element,
    )
    defect = create_slab_defect(configuration=configuration, builder=None)

    assert defect.basis.elements.values[-1] == expected_last_element
    assertion_utils.assert_deep_almost_equal(expected_last_coord, defect.basis.coordinates.values[-1])


@pytest.mark.skip(reason="we'll fix before epic-7623 is merged")
@pytest.mark.parametrize(
    "crystal_config, position_on_surface, distance_z, chemical_element,"
    + " expected_last_element, expected_coords_platform",
    [
        (
            SI_CONVENTIONAL_SLAB_001,
            [0.5, 0.5],
            2,
            "Si",
            "Si",
            {
                "darwin": [5.80049999709975, 3.3489202347599645, 14.234895071861322],
                "other": [5.80049999709975, 3.3489202347599645, 14.234895071861322],
            },
        )
    ],
)
def test_create_adatom_equidistant(
    crystal_config, position_on_surface, distance_z, chemical_element, expected_last_element, expected_coords_platform
):
    # Adatom of Si at approximate 0.5, 0.5 position
    crystal = Material.create(crystal_config)
    configuration = AdatomSlabPointDefectConfiguration(
        crystal=crystal,
        position_on_surface=position_on_surface,
        distance_z=distance_z,
        chemical_element=chemical_element,
    )
    defect = create_slab_defect(configuration=configuration, builder=EquidistantAdatomSlabDefectBuilder())

    assert defect.basis.elements.values[-1] == expected_last_element
    assert (len(crystal.basis.coordinates.values) + 1) == len(defect.basis.coordinates.values)
    defect.to_cartesian()
    # TODO: resolve the problem with the test in GH pipeline
    # on MacOS slab atoms have different coordinates than in GH and pyodide
    # for the same versions of packages
    if sys.platform == "darwin":
        coordinate_expected = expected_coords_platform["darwin"]
    else:
        coordinate_expected = expected_coords_platform["other"]

    defect_coordinate = defect.basis.coordinates.values[-1]
    atol = 10 ** (-COORDINATE_TOLERANCE)
    assertion_utils.assert_deep_almost_equal(coordinate_expected, defect_coordinate, atol=atol)


@pytest.mark.skip(reason="we'll fix before epic-7623 is merged")
@pytest.mark.parametrize(
    "crystal_config, position_on_surface, distance_z, chemical_element,"
    + " expected_last_element, expected_coords_platform",
    [
        (
            SI_CONVENTIONAL_SLAB_001,
            [0.5, 0.5],
            2.5,
            None,
            "Si",
            {"darwin": [0.875, 0.75, 0.534157228], "other": [0.875, 0.75, 0.588188708]},
        )
    ],
)
def test_create_crystal_site_adatom(
    crystal_config, position_on_surface, distance_z, chemical_element, expected_last_element, expected_coords_platform
):
    # Adatom of Si (autodetect) at approximate 0.5, 0.5 position
    crystal = Material.create(crystal_config)
    configuration = AdatomSlabPointDefectConfiguration(
        crystal=crystal,
        position_on_surface=position_on_surface,
        distance_z=distance_z,
        chemical_element=chemical_element,
    )
    builder = CrystalSiteAdatomSlabDefectBuilder()
    defect = create_slab_defect(configuration=configuration, builder=builder)

    assert defect.basis.elements.values[-1] == expected_last_element

    if sys.platform == "darwin":
        coordinates_expected = expected_coords_platform["darwin"]
    else:
        coordinates_expected = expected_coords_platform["other"]
    defect_coordinate = defect.basis.coordinates.values[-1]
    atol = 10 ** (-COORDINATE_TOLERANCE)
    assertion_utils.assert_deep_almost_equal(coordinates_expected, defect_coordinate, atol=atol)


@pytest.mark.parametrize(
    "slab_material_config, position_on_surface, distance_z, expected_center",
    [(SI_CONVENTIONAL_SLAB_001, [0.5, 0.5], 2.5, [0.5, 0.5, 0.7573538315436493])],
)
def test_get_equidistant_position(slab_material_config, position_on_surface, distance_z, expected_center):
    builder = EquidistantAdatomSlabDefectBuilder()

    slab_material = Material.create(slab_material_config)

    equidistant_position = builder.get_equidistant_position(
        material=slab_material, position_on_surface=position_on_surface, distance_z=distance_z
    )

    assertion_utils.assert_deep_almost_equal(equidistant_position, expected_center)
