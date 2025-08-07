from ..base.builder import PointDefectBuilder
from .configuration import (
    InterstitialDefectConfiguration,
)


class InterstitialDefectBuilder(PointDefectBuilder):
    _ConfigurationType = InterstitialDefectConfiguration
