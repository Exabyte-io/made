from mat3ra.made.lattice import COORDINATE_TOLERANCE
from mat3ra.made.material import Material
from mat3ra.made.tools.build.defect import (
    AdatomSlabPointDefectConfiguration,
    CrystalSiteAdatomSlabDefectBuilder,
    EquidistantAdatomSlabDefectBuilder,
    PointDefectBuilderParameters,
    PointDefectConfiguration,
    PointDefectTypeEnum,
    create_defect,
    create_defects,
    create_slab_defect,
)
from mat3ra.made.tools.build.defect.builders import (
    IslandSlabDefectBuilder,
    PointDefectPairBuilder,
    TerraceSlabDefectBuilder,
)
from mat3ra.made.tools.build.defect.configuration import (
    IslandSlabDefectConfiguration,
    PointDefectPairConfiguration,
    TerraceSlabDefectConfiguration,
)
from mat3ra.made.tools.utils import coordinate as CoordinateCondition
from mat3ra.utils import assertion as assertion_utils
from unit.fixtures.generated.fixtures import SLAB_001, SLAB_111

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
        crystal=SLAB_111, position_on_surface=[0.5, 0.5], distance_z=2, chemical_element="Si"
    )
    defect = create_slab_defect(configuration=configuration, builder=None)

    assert defect.basis.elements.values[-1] == "Si"
    assertion_utils.assert_deep_almost_equal([0.5, 0.5, 0.872332562], defect.basis.coordinates.values[-1])


def test_create_adatom_equidistant():
    # Adatom of Si at approximate 0.5, 0.5 position
    configuration = AdatomSlabPointDefectConfiguration(
        crystal=SLAB_111, position_on_surface=[0.5, 0.5], distance_z=2, chemical_element="Si"
    )
    defect = create_slab_defect(configuration=configuration, builder=EquidistantAdatomSlabDefectBuilder())

    assert defect.basis.elements.values[-1] == "Si"
    assert (len(configuration.crystal.basis.coordinates.values) + 1) == len(defect.basis.coordinates.values)
    defect.to_cartesian()
    # TODO: resolve the problem with the test in GH pipeline
    # on MacOS slab atoms have different coordinates than in GH and pyodide
    # for the same versions of packages
    coordinate_macosx = [6.477224996, 3.739627331, 14.234895469]
    coordinate_linux_and_emscripten = [5.123775004, 3.739627331, 14.234895469]
    defect_coordinate = defect.basis.coordinates.values[-1]
    atol = 10 ** (-COORDINATE_TOLERANCE)
    try:
        assertion_utils.assert_deep_almost_equal(coordinate_macosx, defect_coordinate, atol=atol)
    except AssertionError:
        assertion_utils.assert_deep_almost_equal(coordinate_linux_and_emscripten, defect_coordinate, atol=atol)


def test_create_crystal_site_adatom():
    # Adatom of Si (autodetect) at approximate 0.5, 0.5 position
    configuration = AdatomSlabPointDefectConfiguration(
        crystal=SLAB_111, position_on_surface=[0.5, 0.5], distance_z=2, chemical_element=None
    )
    builder = CrystalSiteAdatomSlabDefectBuilder()
    defect = create_slab_defect(configuration=configuration, builder=builder)

    assert defect.basis.elements.values[-1] == "Si"

    coordinates_macosx = [0.458333, 0.458333, 0.628217]
    coordinates_linux_and_emscripten = [0.083333, 0.458333, 0.628217]
    defect_coordinate = defect.basis.coordinates.values[-1]
    atol = 10 ** (-COORDINATE_TOLERANCE)
    try:
        assertion_utils.assert_deep_almost_equal(coordinates_macosx, defect_coordinate, atol=atol)
    except AssertionError:
        assertion_utils.assert_deep_almost_equal(coordinates_linux_and_emscripten, defect_coordinate, atol=atol)


def test_create_island():
    condition = CoordinateCondition.CylinderCoordinateCondition(
        center_position=[0.5, 0.5], radius=0.15, min_z=0, max_z=1
    )
    island_config = IslandSlabDefectConfiguration(
        crystal=SLAB_001,
        defect_type="island",
        condition=condition,
        number_of_added_layers=1,
    )

    defect = create_slab_defect(configuration=island_config, builder=IslandSlabDefectBuilder())

    # Only 1 atoms in the island were added for this configuration with 001 slab orientation
    NUMBER_OF_ATOMS_IN_ISLAND = 1
    assert len(defect.basis.elements.values) == len(SLAB_001.basis.elements.values) + NUMBER_OF_ATOMS_IN_ISLAND
    assert defect.basis.elements.values[-1] == "Si"


def test_create_terrace():
    config = TerraceSlabDefectConfiguration(
        crystal=SLAB_001,
        cut_direction=[1, 0, 0],
        pivot_coordinate=[0.5, 0.5, 0.5],
        number_of_added_layers=1,
    )
    new_slab = TerraceSlabDefectBuilder().get_material(configuration=config)
    coordinate_macosx = [0.777786396, 0.5, 0.414655236]
    coordinate_linux_and_emscripten = [0.627786404, 0.25, 0.439235145]
    defect_coordinate = new_slab.basis.coordinates.values[42]
    atol = 10 ** (-COORDINATE_TOLERANCE)
    try:
        assertion_utils.assert_deep_almost_equal(coordinate_macosx, defect_coordinate, atol=atol)
    except AssertionError:
        assertion_utils.assert_deep_almost_equal(coordinate_linux_and_emscripten, defect_coordinate, atol=atol)


def test_create_defect_pair():
    defect1_config = PointDefectConfiguration.from_approximate_coordinate(
        crystal=SLAB_001,
        defect_type=PointDefectTypeEnum.VACANCY,
        approximate_coordinate=[0.5, 0.5, 0.25],
    )
    defect2_config = PointDefectConfiguration(
        defect_type=PointDefectTypeEnum.INTERSTITIAL, coordinate=[0.5, 0.55, 0.35], chemical_element="P"
    )
    defect_pair_config = PointDefectPairConfiguration(
        primary_defect_configuration=defect1_config,
        secondary_defect_configuration=defect2_config,
    )
    defect_material = PointDefectPairBuilder().get_material(defect_pair_config)

    assertion_utils.assert_deep_almost_equal("P", defect_material.basis.elements.values[0])
    assertion_utils.assert_deep_almost_equal([0.5, 0.55, 0.35], defect_material.basis.coordinates.values[0])
