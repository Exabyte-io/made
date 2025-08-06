from types import SimpleNamespace

import pytest
from mat3ra.made.material import Material
from mat3ra.made.tools.build.defective_structures.zero_dimensional.pair_defect.builder import PairDefectBuilder
from mat3ra.made.tools.build.defective_structures.zero_dimensional.pair_defect.configuration import (
    PairDefectConfiguration,
)
from mat3ra.made.tools.build.defective_structures.zero_dimensional.point_defect import (
    SubstitutionPlacementMethodEnum,
    create_defect_configuration,
)
from unit.fixtures.bulk import BULK_Si_PRIMITIVE
from unit.fixtures.pair_defects import (
    PAIR_DEFECT_SUBSTITUTION_VACANCY_BULK_PRIMITIVE_Si,
    PAIR_DEFECT_VACANCY_INTERSTITIAL_BULK_PRIMITIVE_Si,
)
from unit.utils import assert_two_entities_deep_almost_equal


@pytest.mark.parametrize(
    "material_config, pair_defect_params, expected_material_config",
    [
        (
            BULK_Si_PRIMITIVE,
            SimpleNamespace(
                primary_defect_type="vacancy",
                primary_coordinate=[0.0, 0.0, 0.0],
                secondary_defect_type="interstitial",
                secondary_coordinate=[0.5, 0.5, 0.5],
                secondary_element="Ge",
            ),
            PAIR_DEFECT_VACANCY_INTERSTITIAL_BULK_PRIMITIVE_Si,
        ),
        (
            BULK_Si_PRIMITIVE,
            SimpleNamespace(
                primary_defect_type="substitution",
                primary_coordinate=[0.123, 0.123, 0.123],
                primary_element="Ge",
                primary_placement_method=SubstitutionPlacementMethodEnum.CLOSEST_SITE,
                secondary_defect_type="vacancy",
                secondary_coordinate=[0.25, 0.25, 0.25],
            ),
            PAIR_DEFECT_SUBSTITUTION_VACANCY_BULK_PRIMITIVE_Si,
        ),
    ],
)
def test_pair_defect_builder(material_config, pair_defect_params, expected_material_config):
    crystal = Material.create(material_config)

    primary_config = create_defect_configuration(
        material=crystal,
        defect_type=pair_defect_params.primary_defect_type,
        coordinate=pair_defect_params.primary_coordinate,
        element=pair_defect_params.primary_element if hasattr(pair_defect_params, "primary_element") else None,
        placement_method=(
            pair_defect_params.primary_placement_method
            if hasattr(pair_defect_params, "primary_placement_method")
            else None
        ),
    )
    secondary_config = create_defect_configuration(
        material=crystal,
        defect_type=pair_defect_params.secondary_defect_type,
        coordinate=pair_defect_params.secondary_coordinate,
        element=pair_defect_params.secondary_element if hasattr(pair_defect_params, "secondary_element") else None,
        placement_method=(
            pair_defect_params.secondary_placement_method
            if hasattr(pair_defect_params, "secondary_placement_method")
            else None
        ),
    )

    config = PairDefectConfiguration.from_parameters(
        crystal=crystal,
        primary_defect_configuration=primary_config,
        secondary_defect_configuration=secondary_config,
    )

    builder = PairDefectBuilder()
    defect = builder.get_material(config)

    assert_two_entities_deep_almost_equal(defect, expected_material_config)
