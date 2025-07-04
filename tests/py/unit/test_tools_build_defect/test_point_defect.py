import pytest
from mat3ra.made.material import Material
from mat3ra.made.tools.build.defect import (
    PointDefectBuilderParameters,
    PointDefectConfigurationLegacy,
    PointDefectTypeEnum,
    create_defect,
    create_defects,
)
from mat3ra.made.tools.build.defect.point.helpers import (
    create_vacancy_defect,
    create_substitution_defect,
    create_interstitial_defect,
)
from unit.fixtures.bulk import BULK_Si_PRIMITIVE
from unit.fixtures.point_defects import (
    VACANCY_DEFECT_BULK_PRIMITIVE_Si,
    SUBSTITUTION_DEFECT_BULK_PRIMITIVE_Si,
    INTERSTITIAL_DEFECT_BULK_PRIMITIVE_Si,
)
from unit.utils import assert_two_entities_deep_almost_equal
from types import SimpleNamespace


@pytest.mark.parametrize(
    "crystal_config, defect_type, site_id, expected_elements_len",
    [
        (BULK_Si_PRIMITIVE, "vacancy", 0, 1),
    ],
)
def test_create_vacancy(crystal_config, defect_type, site_id, expected_elements_len):
    # vacancy in place of 0 element
    crystal = Material.create(crystal_config)
    configuration = PointDefectConfigurationLegacy(crystal=crystal, defect_type=defect_type, site_id=site_id)
    defect = create_defect(configuration)

    assert len(defect.basis.elements.values) == expected_elements_len


@pytest.mark.parametrize(
    "crystal_config, defect_type, chemical_element, expected_elements, expected_coordinate",
    [
        (
            BULK_Si_PRIMITIVE,
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
    configuration = PointDefectConfigurationLegacy(
        crystal=crystal, defect_type=defect_type, chemical_element=chemical_element
    )
    defect = create_defect(configuration)

    assert defect.basis.elements.to_dict() == expected_elements
    assert defect.basis.coordinates.to_dict()[0] == expected_coordinate


@pytest.mark.parametrize(
    "crystal_config, defect_type, chemical_element, coordinate, expected_elements",
    [
        (
            BULK_Si_PRIMITIVE,
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
    configuration = PointDefectConfigurationLegacy(
        crystal=crystal, defect_type=defect_type, chemical_element=chemical_element, coordinate=coordinate
    )
    defect = create_defect(configuration)

    assert defect.basis.elements.to_dict() == expected_elements


@pytest.mark.parametrize(
    "crystal_config, defect_type, chemical_element, coordinate, placement_method,"
    + " expected_element, expected_coordinates_platform",
    [
        (
            BULK_Si_PRIMITIVE,
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
    expected_coordinates_platform,
):
    crystal = Material.create(crystal_config)
    configuration = PointDefectConfigurationLegacy(
        crystal=crystal,
        defect_type=defect_type,
        chemical_element=chemical_element,
        # Voronoi must resolve to [0.5, 0.5, 0.5] for Si structure
        coordinate=coordinate,
        placement_method=placement_method,
    )
    defect = create_defect(configuration)
    assert defect.basis.elements.values[-1] == expected_element

    coordinate_x86 = expected_coordinates_platform["x86"]
    coordinate_arm64 = expected_coordinates_platform["arm64"]
    defect_coordinate = defect.basis.coordinates.values[-1]
    is_passing_on_x86 = coordinate_x86 == defect_coordinate
    is_passing_on_arm64 = coordinate_arm64 == defect_coordinate
    assert is_passing_on_x86 or is_passing_on_arm64


@pytest.mark.parametrize(
    "crystal_config, defect_type, chemical_element, site_id, center_defect, expected_elements",
    [
        (
            BULK_Si_PRIMITIVE,
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
    defect_configuration = PointDefectConfigurationLegacy.from_site_id(
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
            BULK_Si_PRIMITIVE,
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
        PointDefectConfigurationLegacy.from_site_id(crystal=crystal, **params) for params in defect_configs_params
    ]
    defect_builder_parameters = PointDefectBuilderParameters(**builder_params) if builder_params else None
    material_with_defect = create_defects(builder_parameters=defect_builder_parameters, configurations=configurations)

    assert material_with_defect.basis.elements.to_dict() == expected_elements


@pytest.mark.parametrize(
    "material_config, defect_params, expected_material_config",
    [
        (
            BULK_Si_PRIMITIVE,
            SimpleNamespace(type="vacancy", coordinate=[0.0, 0.0, 0.0]),
            VACANCY_DEFECT_BULK_PRIMITIVE_Si,
        ),
        (
            BULK_Si_PRIMITIVE,
            SimpleNamespace(type="substitution", coordinate=[0.0, 0.0, 0.0], element="Ge"),
            SUBSTITUTION_DEFECT_BULK_PRIMITIVE_Si,
        ),
        (
            BULK_Si_PRIMITIVE,
            SimpleNamespace(type="interstitial", coordinate=[0.5, 0.5, 0.5], element="C"),
            INTERSTITIAL_DEFECT_BULK_PRIMITIVE_Si,
        ),
    ],
)
def test_point_defect_helpers(material_config, defect_params, expected_material_config):
    from mat3ra.made.tools.build.defect.point.helpers import (
        create_vacancy_defect,
        create_substitution_defect,
        create_interstitial_defect,
    )
    from mat3ra.made.material import Material

    crystal = Material.create(material_config)
    if defect_params.type == "vacancy":
        defect = create_vacancy_defect(crystal, defect_params.coordinate)
    elif defect_params.type == "substitution":
        defect = create_substitution_defect(crystal, defect_params.coordinate, defect_params.element)
    elif defect_params.type == "interstitial":
        defect = create_interstitial_defect(crystal, defect_params.coordinate, defect_params.element)
    else:
        raise ValueError(f"Unknown defect_type: {defect_params.type}")
    assert_two_entities_deep_almost_equal(defect, expected_material_config)
