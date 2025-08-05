from ..base.builder import PointDefectBuilder
from .configuration import SubstitutionalDefectConfiguration


class SubstitutionalDefectBuilder(PointDefectBuilder):
    _ConfigurationType = SubstitutionalDefectConfiguration
