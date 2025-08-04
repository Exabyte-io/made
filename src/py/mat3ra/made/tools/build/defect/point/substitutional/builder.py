from mat3ra.made.tools.build.defect.point.base.builder import PointDefectBuilder
from mat3ra.made.tools.build.defect.point.substitutional.configuration import SubstitutionalDefectConfiguration


class SubstitutionalDefectBuilder(PointDefectBuilder):
    _ConfigurationType = SubstitutionalDefectConfiguration
