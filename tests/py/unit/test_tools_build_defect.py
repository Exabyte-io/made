import pytest
import sys

from mat3ra.made.lattice import COORDINATE_TOLERANCE
from mat3ra.made.material import Material

# fmt: off
from mat3ra.made.tools.build.defect import (
    EquidistantAdatomSlabDefectBuilder,
    PointDefectBuilderParameters,
    PointDefectConfiguration,
    PointDefectTypeEnum,
    create_defect,
    create_defects, create_slab_defect,
)
from mat3ra.made.tools.build.defect.builders import (
    CrystalSiteAdatomSlabDefectBuilder,
    IslandSlabDefectBuilder,
    SlabDefectBuilder,
    SlabDefectBuilderParameters,
    TerraceSlabDefectBuilder,
)
from mat3ra.made.tools.build.defect.configuration import (
    AdatomSlabPointDefectConfiguration,
    IslandSlabDefectConfiguration,
    TerraceSlabDefectConfiguration,
)
from mat3ra.made.tools.utils import coordinate as CoordinateCondition
from mat3ra.utils import assertion as assertion_utils

from unit.fixtures.slab import (
    SI_CONVENTIONAL_SLAB_001,
    SI_SLAB_001_ADDED_FRACTIONAL_LAYER,
    SI_SLAB_001_ADDED_LAYER,
)

# fmt: on


clean_material = Material.create_default()


def test_create_vacancy():
    # vacancy in place of 0 element
    configuration = PointDefectConfiguration(crystal=clean_material, defect_type="vacancy", site_id=0)
    defect = create_defect(configuration)

    assert len(defect.basis.elements.values) == 1


def test_create_substitution():
    # Substitution of Ge in place of Si at default site_id=0
    configuration = PointDefectConfiguration(crystal=clean_material, defect_type="substitution", chemical_element="Ge")
    defect = create_defect(configuration)

    assert defect.basis.elements.to_dict() == [{"id": 0, "value": "Ge"}, {"id": 1, "value": "Si"}]
    assert defect.basis.coordinates.to_dict()[0] == {"id": 0, "value": [0.0, 0.0, 0.0]}


def test_create_interstitial():
    # Interstitial Ge at 0.5, 0.5, 0.5 position
    configuration = PointDefectConfiguration(
        crystal=clean_material, defect_type="interstitial", chemical_element="Ge", coordinate=[0.5, 0.5, 0.5]
    )
    defect = create_defect(configuration)

    assert defect.basis.elements.to_dict() == [
        {"id": 0, "value": "Ge"},
        {"id": 1, "value": "Si"},
        {"id": 2, "value": "Si"},
    ]


def test_create_interstitial_voronoi():
    configuration = PointDefectConfiguration(
        crystal=clean_material,
        defect_type="interstitial",
        chemical_element="Ge",
        # Voronoi must resolve to [0.5, 0.5, 0.5] for Si structure
        coordinate=[0.25, 0.25, 0.5],
        placement_method="voronoi_site",
    )
    defect = create_defect(configuration)
    assert defect.basis.elements.values[-1] == "Ge"

    coordinate_x86 = [0.5, 0.5, 0.5]
    coordinate_arm64 = [0.625, 0.625, 0.125]
    defect_coordinate = defect.basis.coordinates.values[-1]
    is_passing_on_x86 = coordinate_x86 == defect_coordinate
    is_passing_on_arm64 = coordinate_arm64 == defect_coordinate
    assert is_passing_on_x86 or is_passing_on_arm64


def test_create_defect_from_site_id():
    # Substitution of Ge in place of Si at site_id=1
    defect_configuration = PointDefectConfiguration.from_site_id(
        crystal=clean_material, defect_type="substitution", chemical_element="Ge", site_id=1
    )
    defect_builder_parameters = PointDefectBuilderParameters(center_defect=False)
    material_with_defect = create_defect(
        builder_parameters=defect_builder_parameters, configuration=defect_configuration
    )

    assert material_with_defect.basis.elements.to_dict() == [
        {"id": 0, "value": "Si"},
        {"id": 1, "value": "Ge"},
    ]


def test_create_defects():
    # Substitution of Ge in place of Si at site_id=1
    defect_configuration1 = PointDefectConfiguration.from_site_id(
        crystal=clean_material, defect_type=PointDefectTypeEnum.SUBSTITUTION, chemical_element="Ge", site_id=1
    )
    defect_configuration2 = PointDefectConfiguration.from_site_id(
        crystal=clean_material, defect_type=PointDefectTypeEnum.SUBSTITUTION, chemical_element="Ge", site_id=0
    )
    defect_builder_parameters = PointDefectBuilderParameters(center_defect=False)
    material_with_defect = create_defects(
        builder_parameters=defect_builder_parameters, configurations=[defect_configuration1, defect_configuration2]
    )

    assert material_with_defect.basis.elements.to_dict() == [
        {"id": 0, "value": "Ge"},
        {"id": 1, "value": "Ge"},
    ]


def test_create_adatom():
    # Adatom of Si at 0.5, 0.5 position
    configuration = AdatomSlabPointDefectConfiguration(
        crystal=SI_CONVENTIONAL_SLAB_001, position_on_surface=[0.5, 0.5], distance_z=2, chemical_element="Si"
    )
    defect = create_slab_defect(configuration=configuration, builder=None)

    assert defect.basis.elements.values[-1] == "Si"
    assertion_utils.assert_deep_almost_equal([0.5, 0.5, 0.7641022175421057], defect.basis.coordinates.values[-1])


def test_create_adatom_equidistant():
    # Adatom of Si at approximate 0.5, 0.5 position
    configuration = AdatomSlabPointDefectConfiguration(
        crystal=SI_CONVENTIONAL_SLAB_001, position_on_surface=[0.5, 0.5], distance_z=2, chemical_element="Si"
    )
    defect = create_slab_defect(configuration=configuration, builder=EquidistantAdatomSlabDefectBuilder())

    assert defect.basis.elements.values[-1] == "Si"
    assert (len(configuration.crystal.basis.coordinates.values) + 1) == len(defect.basis.coordinates.values)
    defect.to_cartesian()
    # TODO: resolve the problem with the test in GH pipeline
    # on MacOS slab atoms have different coordinates than in GH and pyodide
    # for the same versions of packages
    if sys.platform == "darwin":
        coordinate_expected = [5.80049999709975, 3.3489202347599645, 14.234895071861322]
    else:
        coordinate_expected = [5.80049999709975, 3.3489202347599645, 14.234895071861322]

    defect_coordinate = defect.basis.coordinates.values[-1]
    atol = 10 ** (-COORDINATE_TOLERANCE)
    assertion_utils.assert_deep_almost_equal(coordinate_expected, defect_coordinate, atol=atol)


def test_create_crystal_site_adatom():
    # Adatom of Si (autodetect) at approximate 0.5, 0.5 position
    configuration = AdatomSlabPointDefectConfiguration(
        crystal=SI_CONVENTIONAL_SLAB_001, position_on_surface=[0.5, 0.5], distance_z=2.5, chemical_element=None
    )
    builder = CrystalSiteAdatomSlabDefectBuilder()
    defect = create_slab_defect(configuration=configuration, builder=builder)

    assert defect.basis.elements.values[-1] == "Si"

    if sys.platform == "darwin":
        coordinates_expected = [0.875, 0.75, 0.534157228]
    else:
        coordinates_expected = [0.875, 0.75, 0.588188708]
    defect_coordinate = defect.basis.coordinates.values[-1]
    atol = 10 ** (-COORDINATE_TOLERANCE)
    assertion_utils.assert_deep_almost_equal(coordinates_expected, defect_coordinate, atol=atol)


def test_create_island():
    # TODO: use TiN
    condition = CoordinateCondition.CylinderCoordinateCondition(
        center_position=[0.5, 0.5], radius=0.15, min_z=0, max_z=1
    )
    island_config = IslandSlabDefectConfiguration(
        crystal=SI_CONVENTIONAL_SLAB_001,
        defect_type="island",
        condition=condition,
        number_of_added_layers=1,
    )

    defect = create_slab_defect(configuration=island_config, builder=IslandSlabDefectBuilder())

    # 1 atom in the island were added for this configuration with 001 slab orientation
    NUMBER_OF_ATOMS_IN_ISLAND = 1
    assert (
        len(defect.basis.elements.values)
        == len(island_config.crystal.basis.elements.values) + NUMBER_OF_ATOMS_IN_ISLAND
    )
    assert defect.basis.elements.values[-1] == "Si"


def test_create_terrace():
    config = TerraceSlabDefectConfiguration(
        crystal=SI_CONVENTIONAL_SLAB_001,
        cut_direction=[1, 0, 0],
        pivot_coordinate=[0.5, 0.5, 0.5],
        number_of_added_layers=1,
    )
    new_slab = TerraceSlabDefectBuilder().get_material(configuration=config)
    if sys.platform == "darwin":
        coordinate_expected = [0.627786405, 0.75, 0.671264194]
    else:
        coordinate_expected = [0.591583068, 0.75, 0.716534426]
    defect_coordinate = new_slab.basis.coordinates.values[-1]  # Use last atom (index 59) instead of previous index 35
    atol = 10 ** (-COORDINATE_TOLERANCE)
    assertion_utils.assert_deep_almost_equal(coordinate_expected, defect_coordinate, atol=atol)


@pytest.mark.skip(reason="Slab with additional layers should be adjusted")
def test_create_material_with_additional_layers():
    """Test adding layers to a slab material"""
    # Create the builder
    builder_params = SlabDefectBuilderParameters(auto_add_vacuum=True, vacuum_thickness=5.0)
    builder = SlabDefectBuilder(build_parameters=builder_params)

    # Test adding 1 layer to SI_SLAB_001
    original_slab = Material.create(SI_CONVENTIONAL_SLAB_001)
    slab_with_additional_layer = builder.create_material_with_additional_layers(original_slab, 1)

    assertion_utils.assert_deep_almost_equal(slab_with_additional_layer, SI_SLAB_001_ADDED_LAYER)


@pytest.mark.skip(reason="Slab with additional layers should be adjusted")
def test_create_material_with_additional_fractional_layers():
    """Test adding fractional layers to a slab material"""
    # Create the builder
    builder_params = SlabDefectBuilderParameters(auto_add_vacuum=True, vacuum_thickness=5.0)
    builder = SlabDefectBuilder(build_parameters=builder_params)

    # Test adding 1.5 layers to SI_SLAB_001
    original_slab = Material.create(SI_CONVENTIONAL_SLAB_001)
    slab_with_fractional_layer = builder.create_material_with_additional_layers(original_slab, 1.5)

    # Compare with expected fixture
    assertion_utils.assert_deep_almost_equal(slab_with_fractional_layer, SI_SLAB_001_ADDED_FRACTIONAL_LAYER)


def test_get_equidistant_position():
    builder = EquidistantAdatomSlabDefectBuilder()

    slab_material = Material.create(SI_CONVENTIONAL_SLAB_001)

    position_on_surface = [0.5, 0.5]
    distance_z = 2.5

    equidistant_position = builder.get_equidistant_position(
        material=slab_material, position_on_surface=position_on_surface, distance_z=distance_z
    )

    expected_center = [0.5, 0.5, 0.7573538315436493]
    assertion_utils.assert_deep_almost_equal(equidistant_position, expected_center)
