from mat3ra.made.tools.build.defect.point.base.builder import PointDefectBuilder
from mat3ra.made.tools.build.defect.point.intersitital.configuration import (
    InterstitialDefectConfiguration,
)


class InterstitialDefectBuilder(PointDefectBuilder):
    _ConfigurationType = InterstitialDefectConfiguration
