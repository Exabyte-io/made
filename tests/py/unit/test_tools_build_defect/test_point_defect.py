from types import SimpleNamespace

import pytest

from mat3ra.made.material import Material
from mat3ra.made.tools.build.defect.point.helpers import (
    create_interstitial_defect,
    create_substitution_defect,
    create_vacancy_defect,
)
from mat3ra.made.tools.build.defect.enums import AtomPlacementMethodEnum
from unit.fixtures.bulk import BULK_Si_PRIMITIVE
from unit.fixtures.point_defects import (
    INTERSTITIAL_DEFECT_BULK_PRIMITIVE_Si,
    SUBSTITUTION_DEFECT_BULK_PRIMITIVE_Si,
    VACANCY_DEFECT_BULK_PRIMITIVE_Si,
    INTERSTITIAL_VORONOI_DEFECT_BULK_PRIMITIVE_Si,
)
from unit.utils import assert_two_entities_deep_almost_equal


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
        (
            BULK_Si_PRIMITIVE,
            SimpleNamespace(
                type="interstitial",
                coordinate=[0.25, 0.25, 0.5],
                element="Ge",
                placement_method=AtomPlacementMethodEnum.VORONOI_SITE,
            ),
            INTERSTITIAL_VORONOI_DEFECT_BULK_PRIMITIVE_Si,
        ),
    ],
)
def test_point_defect_helpers(material_config, defect_params, expected_material_config):
    crystal = Material.create(material_config)

    if defect_params.type == "vacancy":
        placement_method = getattr(defect_params, "placement_method", AtomPlacementMethodEnum.COORDINATE)
        defect = create_vacancy_defect(crystal, defect_params.coordinate, placement_method)
    elif defect_params.type == "substitution":
        placement_method = getattr(defect_params, "placement_method", AtomPlacementMethodEnum.COORDINATE)
        defect = create_substitution_defect(crystal, defect_params.coordinate, defect_params.element, placement_method)
    elif defect_params.type == "interstitial":
        placement_method = getattr(defect_params, "placement_method", AtomPlacementMethodEnum.COORDINATE)
        defect = create_interstitial_defect(crystal, defect_params.coordinate, defect_params.element, placement_method)
    else:
        raise ValueError(f"Unknown defect_type: {defect_params.type}")

    assert_two_entities_deep_almost_equal(defect, expected_material_config)
