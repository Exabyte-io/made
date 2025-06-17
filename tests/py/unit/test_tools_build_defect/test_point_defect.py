import sys

import pytest
from mat3ra.made.material import Material
from mat3ra.made.tools.build.defect import (
    PointDefectBuilderParameters,
    PointDefectConfiguration,
    PointDefectTypeEnum,
    create_defect,
    create_defects,
)


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
    "crystal_config, defect_type, chemical_element, coordinate, placement_method,"
    + " expected_element, expected_coords_platform",
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