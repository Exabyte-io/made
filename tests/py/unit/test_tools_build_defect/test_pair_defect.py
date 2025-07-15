from mat3ra.made.material import Material
from mat3ra.made.tools.build.defect.pair_defect.builders import PairDefectBuilder
from mat3ra.made.tools.build.defect.pair_defect.configuration import PairDefectConfiguration

from mat3ra.made.tools.build.defect.point.configuration import (
    VacancyDefectConfiguration,
    SubstitutionalDefectConfiguration,
)
from mat3ra.made.tools.build.defect.point.helpers import create_multiple_defects
from unit.fixtures.bulk import BULK_Si_PRIMITIVE


def test_pair_defect_builder():
    crystal = Material.create(BULK_Si_PRIMITIVE)

    primary_config = VacancyDefectConfiguration.from_parameters(crystal=crystal, coordinate=[0.0, 0.0, 0.0])
    secondary_config = SubstitutionalDefectConfiguration.from_parameters(
        crystal=crystal, coordinate=[0.5, 0.5, 0.5], element="Ge"
    )

    config = PairDefectConfiguration.from_parameters(
        crystal=crystal,
        primary_defect_configuration=primary_config,
        secondary_defect_configuration=secondary_config,
    )

    builder = PairDefectBuilder()
    defect = builder.get_material(config)

    assert defect is not None


def test_create_multiple_defects():
    material = Material.create(BULK_Si_PRIMITIVE)

    defect_configs = [
        VacancyDefectConfiguration.from_parameters(crystal=material, coordinate=[0.0, 0.0, 0.0]),
        SubstitutionalDefectConfiguration.from_parameters(crystal=material, coordinate=[0.5, 0.5, 0.5], element="Ge"),
    ]

    result = create_multiple_defects(
        material=material,
        defect_configurations=defect_configs,
    )

    assert result is not None
