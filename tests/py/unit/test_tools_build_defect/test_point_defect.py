from types import SimpleNamespace

import pytest
from mat3ra.made.material import Material
from mat3ra.made.tools.build import MaterialWithBuildMetadata
from mat3ra.made.tools.build.defective_structures.zero_dimensional.point_defect import (
    InterstitialPlacementMethodEnum,
    SubstitutionPlacementMethodEnum,
    VacancyPlacementMethodEnum,
)
from mat3ra.made.tools.build.defective_structures.zero_dimensional.point_defect.helpers import (
    create_defect_point_interstitial,
    create_defect_point_substitution,
    create_defect_point_vacancy,
    create_multiple_defects,
)
from mat3ra.made.tools.build.defective_structures.zero_dimensional.point_defect.types import PointDefectDict
from unit.fixtures.bulk import BULK_Si_CONVENTIONAL, BULK_Si_PRIMITIVE
from unit.fixtures.point_defects import (
    INTERSTITIAL_DEFECT_BULK_PRIMITIVE_Si,
    INTERSTITIAL_VORONOI_DEFECT_BULK_PRIMITIVE_Si,
    MULTIPLE_POINT_DEFECTS_BULK_Si_CONVENTIONAL,
    SUBSTITUTION_DEFECT_BULK_PRIMITIVE_Si,
    VACANCY_DEFECT_BULK_PRIMITIVE_Si,
)
from unit.utils import assert_two_entities_deep_almost_equal, get_platform_specific_value


@pytest.mark.parametrize(
    "material_config, defect_params, expected_material_config",
    [
        (
            BULK_Si_PRIMITIVE,
            SimpleNamespace(
                type="vacancy",
                coordinate=[0.0, 0.0, 0.0],
                placement_method=VacancyPlacementMethodEnum.CLOSEST_SITE.value,
            ),
            VACANCY_DEFECT_BULK_PRIMITIVE_Si,
        ),
        (
            BULK_Si_PRIMITIVE,
            SimpleNamespace(
                type="substitution",
                coordinate=[0.0, 0.0, 0.0],
                element="Ge",
                placement_method=SubstitutionPlacementMethodEnum.CLOSEST_SITE.value,
            ),
            SUBSTITUTION_DEFECT_BULK_PRIMITIVE_Si,
        ),
        (
            BULK_Si_PRIMITIVE,
            SimpleNamespace(
                type="interstitial",
                coordinate=[0.5, 0.5, 0.5],
                element="C",
                placement_method=InterstitialPlacementMethodEnum.EXACT_COORDINATE.value,
            ),
            INTERSTITIAL_DEFECT_BULK_PRIMITIVE_Si,
        ),
        (
            BULK_Si_PRIMITIVE,
            SimpleNamespace(
                type="interstitial",
                coordinate=[0.25, 0.25, 0.5],
                element="Ge",
                placement_method=InterstitialPlacementMethodEnum.VORONOI_SITE.value,
            ),
            INTERSTITIAL_VORONOI_DEFECT_BULK_PRIMITIVE_Si,
        ),
    ],
)
def test_point_defect_helpers(material_config, defect_params, expected_material_config):
    crystal = Material.create(material_config)

    if defect_params.type == "vacancy":
        defect = create_defect_point_vacancy(crystal, defect_params.coordinate, defect_params.placement_method)
    elif defect_params.type == "substitution":
        defect = create_defect_point_substitution(
            crystal, defect_params.coordinate, defect_params.element, defect_params.placement_method
        )
    elif defect_params.type == "interstitial":
        defect = create_defect_point_interstitial(
            crystal, defect_params.coordinate, defect_params.element, defect_params.placement_method
        )
    else:
        raise ValueError(f"Unknown defect_type: {defect_params.type}")

    expected_material_config = get_platform_specific_value(expected_material_config)

    assert_two_entities_deep_almost_equal(defect, expected_material_config)


@pytest.mark.parametrize(
    "material_config, defect_params_list, expected_material_config",
    [
        (
            BULK_Si_CONVENTIONAL,
            [
                SimpleNamespace(
                    defect_type="vacancy",
                    coordinate=[0.75, 0.70, 0.70],
                    placement_method=VacancyPlacementMethodEnum.CLOSEST_SITE.value,
                ),
                SimpleNamespace(
                    defect_type="interstitial",
                    coordinate=[0.25, 0.25, 0.2],
                    element="N",
                    placement_method=InterstitialPlacementMethodEnum.VORONOI_SITE.value,
                ),
                SimpleNamespace(
                    defect_type="substitution",
                    coordinate=[0.543, 0.543, 0.5],
                    element="Ge",
                    placement_method="closest_site",
                ),
            ],
            MULTIPLE_POINT_DEFECTS_BULK_Si_CONVENTIONAL,
        ),
    ],
)
def test_create_multiple_defects(material_config, defect_params_list, expected_material_config):
    material = MaterialWithBuildMetadata.create(material_config)

    defect_dicts = []
    for defect_params in defect_params_list:
        defect_data = {
            "type": defect_params.defect_type,
            "coordinate": defect_params.coordinate,
            "placement_method": defect_params.placement_method if hasattr(defect_params, "placement_method") else None,
            "element": defect_params.element if hasattr(defect_params, "element") else None,
        }

        defect_dict = PointDefectDict(**defect_data)

        defect_dicts.append(defect_dict)

    defects = create_multiple_defects(
        material=material,
        defect_dicts=defect_dicts,
    )
    assert_two_entities_deep_almost_equal(defects, expected_material_config)
