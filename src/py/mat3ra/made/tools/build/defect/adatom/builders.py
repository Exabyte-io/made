from .configuration import AdatomDefectConfiguration
from ..slab.builders import SlabStackBuilder


class AdatomDefectBuilder(SlabStackBuilder):
    _ConfigurationType = AdatomDefectConfiguration
