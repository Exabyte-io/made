from types import SimpleNamespace

import pytest

from mat3ra.made.material import Material
from mat3ra.made.tools.build.defect import (
    PointDefectConfigurationLegacy,
    create_defect,
)
from mat3ra.made.tools.build.defect.point.helpers import (
    create_interstitial_defect,
    create_substitution_defect,
    create_vacancy_defect,
)
from unit.fixtures.bulk import BULK_Si_PRIMITIVE
from unit.fixtures.point_defects import (
    INTERSTITIAL_DEFECT_BULK_PRIMITIVE_Si,
    SUBSTITUTION_DEFECT_BULK_PRIMITIVE_Si,
    VACANCY_DEFECT_BULK_PRIMITIVE_Si,
)
from unit.utils import assert_two_entities_deep_almost_equal


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
    defect = create_defect(
    assert defect.basis.elements.values[-1] == expected_element

    coordinate_x86 = expected_coordinates_platform["x86"]
    coordinate_arm64 = expected_coordinates_platform["arm64"]
    defect_coordinate = defect.basis.coordinates.values[-1]
    is_passing_on_x86 = coordinate_x86 == defect_coordinate
    is_passing_on_arm64 = coordinate_arm64 == defect_coordinate
    assert is_passing_on_x86 or is_passing_on_arm64


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
