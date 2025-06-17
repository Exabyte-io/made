import sys

import pytest
from mat3ra.utils import assertion as assertion_utils

from mat3ra.made.lattice import COORDINATE_TOLERANCE
from mat3ra.made.material import Material
from mat3ra.made.tools.build.defect import (
    EquidistantAdatomSlabDefectBuilder,
    PointDefectBuilderParameters,
    PointDefectConfiguration,
    PointDefectTypeEnum,
    create_defect,
    create_defects,
    create_slab_defect,
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
from unit.fixtures.slab import SI_CONVENTIONAL_SLAB_001, SI_SLAB_001_ADDED_FRACTIONAL_LAYER, SI_SLAB_001_ADDED_LAYER

clean_material = Material.create_default()


@pytest.mark.parametrize(
    "crystal_config, defect_type, site_id, expected_elements_len",
    [
        (Material.__default_config__, "vacancy", 0, 1),
    ],
)
def test_create_vacancy(crystal_config, defect_type, site_id, expected_elements_len):
    # vacancy in place of 0 element
    crystal = Material.create(crystal_config)
    configuration = PointDefectConfiguration(crystal=crystal, defect_type=defect_type, site_id=site_id)
    defect = create_defect(configuration)

    assert len(defect.basis.elements.values) == expected_elements_len


@pytest.mark.parametrize(
    "crystal_config, defect_type, chemical_element, expected_elements, expected_coordinate",
    [
        (
            Material.__default_config__,
            "substitution",
            "Ge",
            [{"id": 0, "value": "Ge"}, {"id": 1, "value": "Si"}],
            {"id": 0, "value": [0.0, 0.0, 0.0]},
        ),
    ],
)
def test_create_substitution(crystal_config, defect_type, chemical_element, expected_elements, expected_coordinate):
    # Substitution of Ge in place of Si at default site_id=0
    crystal = Material.create(crystal_config)
    configuration = PointDefectConfiguration(
        crystal=crystal, defect_type=defect_type, chemical_element=chemical_element
    )
    defect = create_defect(configuration)

    assert defect.basis.elements.to_dict() == expected_elements
    assert defect.basis.coordinates.to_dict()[0] == expected_coordinate


@pytest.mark.parametrize(
    "crystal_config, defect_type, chemical_element, coordinate, expected_elements",
    [
        (
            Material.__default_config__,
            "interstitial",
            "Ge",
            [0.5, 0.5, 0.5],
            [{"id": 0, "value": "Ge"}, {"id": 1, "value": "Si"}, {"id": 2, "value": "Si"}],
        ),
    ],
)
def test_create_interstitial(crystal_config, defect_type, chemical_element, coordinate, expected_elements):
    # Interstitial Ge at 0.5, 0.5, 0.5 position
    crystal = Material.create(crystal_config)
    configuration = PointDefectConfiguration(
        crystal=crystal, defect_type=defect_type, chemical_element=chemical_element, coordinate=coordinate
    )
    defect = create_defect(configuration)

    assert defect.basis.elements.to_dict() == expected_elements


@pytest.mark.parametrize(
    "crystal_config, defect_type, chemical_element, coordinate, placement_method, expected_element, expected_coords_platform",
    [
        (
            Material.__default_config__,
            "interstitial",
            "Ge",
            [0.25, 0.25, 0.5],
            "voronoi_site",
            "Ge",
            {"x86": [0.5, 0.5, 0.5], "arm64": [0.625, 0.625, 0.125]},
        ),
    ],
)
def test_create_interstitial_voronoi(
    crystal_config,
    defect_type,
    chemical_element,
    coordinate,
    placement_method,
    expected_element,
    expected_coords_platform,
):
    crystal = Material.create(crystal_config)
    configuration = PointDefectConfiguration(
        crystal=crystal,
        defect_type=defect_type,
        chemical_element=chemical_element,
        # Voronoi must resolve to [0.5, 0.5, 0.5] for Si structure
        coordinate=coordinate,
        placement_method=placement_method,
    )
    defect = create_defect(configuration)
    assert defect.basis.elements.values[-1] == expected_element

    coordinate_x86 = expected_coords_platform["x86"]
    coordinate_arm64 = expected_coords_platform["arm64"]
    defect_coordinate = defect.basis.coordinates.values[-1]
    is_passing_on_x86 = coordinate_x86 == defect_coordinate
    is_passing_on_arm64 = coordinate_arm64 == defect_coordinate
    assert is_passing_on_x86 or is_passing_on_arm64


@pytest.mark.parametrize(
    "crystal_config, defect_type, chemical_element, site_id, center_defect, expected_elements",
    [
        (
            Material.__default_config__,
            "substitution",
            "Ge",
            1,
            False,
            [{"id": 0, "value": "Si"}, {"id": 1, "value": "Ge"}],
        )
    ],
)
def test_create_defect_from_site_id(
    crystal_config, defect_type, chemical_element, site_id, center_defect, expected_elements
):
    # Substitution of Ge in place of Si at site_id=1
    crystal = Material.create(crystal_config)
    defect_configuration = PointDefectConfiguration.from_site_id(
        crystal=crystal, defect_type=defect_type, chemical_element=chemical_element, site_id=site_id
    )
    defect_builder_parameters = PointDefectBuilderParameters(center_defect=center_defect)
    material_with_defect = create_defect(
        builder_parameters=defect_builder_parameters, configuration=defect_configuration
    )

    assert material_with_defect.basis.elements.to_dict() == expected_elements


@pytest.mark.parametrize(
    "crystal_config, defect_configs_params, builder_params, expected_elements",
    [
        (
            Material.__default_config__,
            [
                {"defect_type": PointDefectTypeEnum.SUBSTITUTION, "chemical_element": "Ge", "site_id": 1},
                {"defect_type": PointDefectTypeEnum.SUBSTITUTION, "chemical_element": "Ge", "site_id": 0},
            ],
            {"center_defect": False},
            [{"id": 0, "value": "Ge"}, {"id": 1, "value": "Ge"}],
        )
    ],
)
def test_create_defects(crystal_config, defect_configs_params, builder_params, expected_elements):
    # Substitution of Ge in place of Si at site_id=1
    crystal = Material.create(crystal_config)
    configurations = [
        PointDefectConfiguration.from_site_id(crystal=crystal, **params) for params in defect_configs_params
    ]
    defect_builder_parameters = PointDefectBuilderParameters(**builder_params) if builder_params else None
    material_with_defect = create_defects(builder_parameters=defect_builder_parameters, configurations=configurations)

    assert material_with_defect.basis.elements.to_dict() == expected_elements


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
    "crystal_config, position_on_surface, distance_z, chemical_element, expected_last_element, expected_coords_platform",
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
    "crystal_config, position_on_surface, distance_z, chemical_element, expected_last_element, expected_coords_platform",
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


@pytest.mark.skip(reason="we'll fix before epic-7623 is merged")
@pytest.mark.parametrize(
    "crystal_config, condition_params, defect_type, num_added_layers, num_atoms_in_island, expected_last_element",
    [
        (
            SI_CONVENTIONAL_SLAB_001,
            {"center_position": [0.5, 0.5], "radius": 0.15, "min_z": 0, "max_z": 1},
            "island",
            1,
            1,
            "Si",
        )
    ],
)
def test_create_island(
    crystal_config, condition_params, defect_type, num_added_layers, num_atoms_in_island, expected_last_element
):
    # TODO: use TiN
    crystal = Material.create(crystal_config)
    condition = CoordinateCondition.CylinderCoordinateCondition(**condition_params)
    island_config = IslandSlabDefectConfiguration(
        crystal=crystal,
        defect_type=defect_type,
        condition=condition,
        number_of_added_layers=num_added_layers,
    )

    defect = create_slab_defect(configuration=island_config, builder=IslandSlabDefectBuilder())

    # 1 atom in the island were added for this configuration with 001 slab orientation
    assert len(defect.basis.elements.values) == len(crystal.basis.elements.values) + num_atoms_in_island
    assert defect.basis.elements.values[-1] == expected_last_element


@pytest.mark.skip(reason="we'll fix before epic-7623 is merged")
@pytest.mark.parametrize(
    "crystal_config, cut_direction, pivot_coordinate, num_added_layers, expected_coords_platform",
    [
        (
            SI_CONVENTIONAL_SLAB_001,
            [1, 0, 0],
            [0.5, 0.5, 0.5],
            1,
            {"darwin": [0.627786405, 0.75, 0.671264194], "other": [0.591583068, 0.75, 0.716534426]},
        )
    ],
)
def test_create_terrace(crystal_config, cut_direction, pivot_coordinate, num_added_layers, expected_coords_platform):
    crystal = Material.create(crystal_config)
    config = TerraceSlabDefectConfiguration(
        crystal=crystal,
        cut_direction=cut_direction,
        pivot_coordinate=pivot_coordinate,
        number_of_added_layers=num_added_layers,
    )
    new_slab = TerraceSlabDefectBuilder().get_material(configuration=config)
    if sys.platform == "darwin":
        coordinate_expected = expected_coords_platform["darwin"]
    else:
        coordinate_expected = expected_coords_platform["other"]
    defect_coordinate = new_slab.basis.coordinates.values[-1]  # Use last atom (index 59) instead of previous index 35
    atol = 10 ** (-COORDINATE_TOLERANCE)
    assertion_utils.assert_deep_almost_equal(coordinate_expected, defect_coordinate, atol=atol)


@pytest.mark.skip(reason="Slab with additional layers should be adjusted")
@pytest.mark.parametrize(
    "original_slab_config, layers_to_add, builder_params_dict, expected_slab_config",
    [(SI_CONVENTIONAL_SLAB_001, 1, {"auto_add_vacuum": True, "vacuum_thickness": 5.0}, SI_SLAB_001_ADDED_LAYER)],
)
def test_create_material_with_additional_layers(
    original_slab_config, layers_to_add, builder_params_dict, expected_slab_config
):
    """Test adding layers to a slab material"""
    # Create the builder
    builder_params = SlabDefectBuilderParameters(**builder_params_dict)
    builder = SlabDefectBuilder(build_parameters=builder_params)

    # Test adding 1 layer to SI_SLAB_001
    original_slab = Material.create(original_slab_config)
    expected_slab = Material.create(expected_slab_config)
    slab_with_additional_layer = builder.create_material_with_additional_layers(original_slab, layers_to_add)

    assertion_utils.assert_deep_almost_equal(slab_with_additional_layer, expected_slab)


@pytest.mark.skip(reason="Slab with additional layers should be adjusted")
@pytest.mark.parametrize(
    "original_slab_config, layers_to_add, builder_params_dict, expected_slab_config",
    [
        (
            SI_CONVENTIONAL_SLAB_001,
            1.5,
            {"auto_add_vacuum": True, "vacuum_thickness": 5.0},
            SI_SLAB_001_ADDED_FRACTIONAL_LAYER,
        )
    ],
)
def test_create_material_with_additional_fractional_layers(
    original_slab_config, layers_to_add, builder_params_dict, expected_slab_config
):
    """Test adding fractional layers to a slab material"""
    # Create the builder
    builder_params = SlabDefectBuilderParameters(**builder_params_dict)
    builder = SlabDefectBuilder(build_parameters=builder_params)

    # Test adding 1.5 layers to SI_SLAB_001
    original_slab = Material.create(original_slab_config)
    expected_slab = Material.create(expected_slab_config)
    slab_with_fractional_layer = builder.create_material_with_additional_layers(original_slab, layers_to_add)

    # Compare with expected fixture
    assertion_utils.assert_deep_almost_equal(slab_with_fractional_layer, expected_slab)


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
